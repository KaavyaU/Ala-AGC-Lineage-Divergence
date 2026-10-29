"""
tRNA Ala-AGC Phylogenetic Tree Visualizer
==========================================
Builds a phylogenetic tree from NCBI taxids, colored by A54/U54 lineage.

Requirements:
    pip install ete3 pandas requests

Usage:
    python trna_phylo_tree.py --csv your_data.csv --output tree.svg
    python trna_phylo_tree.py --csv your_data.csv --output tree.png

CSV columns expected:
    species | ncbi taxid | total ala-agc count | A54 count | U54 count | lineage

Lineage values:
    "A54"  -> colored red/coral
    "U54"  -> colored blue
    "Both" -> pie chart showing A54:U54 ratio
    "N/A"  -> excluded
"""

import argparse
import ssl
import sys

import pandas as pd
from ete3 import (
    Tree, TreeStyle, NodeStyle,
    faces, CircleFace, TextFace,
)


def configure_ncbi_ssl() -> None:
    """Work around self-signed certificates in some local Python environments."""
    try:
        ssl._create_default_https_context = ssl._create_unverified_context
    except AttributeError:
        pass

# --- Colors ------------------------------------------------------
COLOR_A54   = "#E8593C"   # warm coral/red  → A54
COLOR_U54   = "#3B8BD4"   # clear blue      → U54
COLOR_OTHER  = "#F0B429"   # gold/amber → atypical (C, G, etc.)
COLOR_LABEL = "#2C2C2A"
COLOR_INNER = "#888780"   # inner node dots

# --- NCBI taxonomy helpers ------------------------------------------------------

def normalize_column_name(name: str) -> str:
    """Collapse punctuation and whitespace so column aliases are easier to match."""
    return "".join(ch.lower() for ch in str(name) if ch.isalnum())


def build_tree_from_taxids(taxids: list[int]) -> Tree:
    """
    Use ETE3's NCBITaxa to build a topology from a list of taxids.
    Returns a Tree object with taxid as node names.
    """
    from ete3 import NCBITaxa
    ncbi = NCBITaxa()

    print("Building tree topology from NCBI taxonomy...")
    try:
        tree = ncbi.get_topology(taxids, intermediate_nodes=True)
        return tree
    except Exception as e:
        print(f"Error building NCBI topology: {e}")
        print("Tip: try running `ncbi.update_taxonomy_database()` once to download the db.")
        sys.exit(1)


# --- Node face: pie chart ------------------------------------------------------

def make_leaf_face(row: pd.Series, size: int = 30) -> faces.Face:
    """
    Returns a StackedBarFace for all lineages (A54, U54, Both),
    proportionally split into A54 | U54 | Other segments.
    """
    a54   = max(row["A54 count"], 0)
    u54   = max(row["U54 count"], 0)
    total = max(row.get("total ala-agc count", a54 + u54), 0)
    other = max(total - a54 - u54, 0)
 
    if total == 0:
        return CircleFace(radius=size // 4, color="#bbbbbb", style="circle")
 
    pct_a54   = 100 * a54   / total
    pct_u54   = 100 * u54   / total
    pct_other = 100 * other / total
 
    bar = faces.StackedBarFace(
        percents=[pct_a54, pct_u54, pct_other],
        width=size * 3,
        height=size // 2,
        colors=[COLOR_A54, COLOR_U54, COLOR_OTHER],
        line_color="#444444",
    )
    bar.margin_left  = 4
    bar.margin_right = 6
    return bar


# --- Tree styling ------------------------------------------------------

def style_tree(tree: Tree, df: pd.DataFrame) -> TreeStyle:
    """Apply node styles and faces; return a configured TreeStyle."""
 
    # Index df by taxid for fast lookup
    df_indexed = df.set_index("ncbi taxid")
 
    from ete3 import NCBITaxa
    ncbi = NCBITaxa()
 
    # Map taxid → scientific name for all nodes
    all_taxids = [int(n.name) for n in tree.traverse() if n.name.isdigit()]
    taxid2name = {}
    try:
        taxid2name = ncbi.get_taxid_translator(all_taxids)
    except Exception:
        pass
 
    for node in tree.traverse():
        node.img_style["size"] = 0  # hide default dot
 
        taxid = int(node.name) if node.name.isdigit() else None
 
        if node.is_leaf() and taxid and taxid in df_indexed.index:
            row = df_indexed.loc[taxid]
 
            # Leaf symbol (stacked bar or circle)
            face = make_leaf_face(row)
            node.add_face(face, column=0, position="branch-right")
 
            # Species name label
            species_name = str(row.get("species", taxid2name.get(taxid, str(taxid))))
            name_face = TextFace(
                f"  {species_name}",
                fsize=9,
                fgcolor=COLOR_LABEL,
                bold=False,
            )
            node.add_face(name_face, column=1, position="branch-right")
 
            # Count annotation (small, muted)
            lineage = str(row["lineage"]).strip()
            if lineage == "Both":
                a54   = int(row["A54 count"])
                u54   = int(row["U54 count"])
                total = int(row.get("total ala-agc count", a54 + u54))
                other = max(total - a54 - u54, 0)
                count_text = f"  A54:{a54}  U54:{u54}" + (f"  Other:{other}" if other else "")
            elif lineage == "A54":
                count_text = f"  A54:{int(row['A54 count'])}"
            elif lineage == "U54":
                count_text = f"  U54:{int(row['U54 count'])}"
            else:
                count_text = ""
 
            if count_text:
                count_face = TextFace(count_text, fsize=7, fgcolor="#888780")
                node.add_face(count_face, column=2, position="branch-right")
 
        else:
            # Internal node: small gray dot
            ns = NodeStyle()
            ns["size"] = 4
            ns["fgcolor"] = COLOR_INNER
            ns["hz_line_width"] = 1
            ns["vt_line_width"] = 1
            node.set_style(ns)
 
            # Optionally label internal nodes with rank/name
            if taxid and taxid in taxid2name:
                sci = taxid2name[taxid]
                inner_face = TextFace(f" {sci}", fsize=7, fgcolor="#aaaaaa")
                node.add_face(inner_face, column=0, position="branch-top")

    # --- Legend ------------------------------------------------------ 
    ts = TreeStyle()
    ts.show_leaf_name = False       # we add our own name faces
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.mode = "r"                   # rectangular (r) or circular (c)
    ts.scale = 40                   # pixels per branch length unit
    ts.branch_vertical_margin = 6
    ts.margin_top    = 30
    ts.margin_bottom = 30
    ts.margin_left   = 20
    ts.margin_right  = 60
 
    # Title
    title = TextFace("tRNA Ala-AGC · Position 54 Lineage", fsize=14,
                     fgcolor=COLOR_LABEL, bold=True)
    title.margin_bottom = 10
    ts.title.add_face(title, column=0)
 
    # Legend entries
    legend_items = [
        (COLOR_A54,   "A54 lineage"),
        (COLOR_U54,   "U54 lineage"),
        (COLOR_OTHER, "Atypical (C / G / other)"),
    ]
    for i, (color, label) in enumerate(legend_items):
        cf = CircleFace(8, color)
        lf = TextFace(f"  {label}", fsize=8, fgcolor=COLOR_LABEL)
        ts.legend.add_face(cf, column=0)
        ts.legend.add_face(lf, column=1)
 
    ts.legend_position = 1  # top-left
 
    return ts


# --- Main ------------------------------------------------------

def main():
    configure_ncbi_ssl()

    parser = argparse.ArgumentParser(description="Build phylogenetic tree from tRNA lineage CSV")
    parser.add_argument("--csv",    required=True,  help="Path to input CSV file")
    parser.add_argument("--output", required=True,  help="Output path (.svg or .png)")
    parser.add_argument("--mode",   default="r",    choices=["r", "c"],
                        help="Tree mode: r=rectangular (default), c=circular")
    parser.add_argument("--update-db", action="store_true",
                        help="Re-download NCBI taxonomy database before running")
    args = parser.parse_args()
 
    # --- Load & clean CSV ------------------------------------------------------
    print(f"Loading {args.csv}...")
    df = pd.read_csv(args.csv)

    normalized_columns = {str(col): normalize_column_name(col) for col in df.columns}
    df = df.rename(columns={col: normalized_columns[col] for col in df.columns})

    col_map = {
        "ncbi taxid": ["ncbitaxid", "taxid", "ncbitaxid"],
        "species": ["species", "organism", "speciesname"],
        "total ala-agc count": [
            "totalalaagccount",
            "totalcount",
            "totalalaagccopies",
            "totalalaagccopiesscore50",
            "totalalaagccopiesscore",
        ],
        "A54 count": ["a54count", "a54", "a54copies", "a54copiesscore50"],
        "U54 count": ["u54count", "u54", "u54copies", "u54copiesscore50"],
        "lineage": ["lineage", "type", "category"],
    }
    for canonical, aliases in col_map.items():
        for alias in aliases:
            if alias in df.columns and canonical.lower() != alias:
                df = df.rename(columns={alias: canonical.lower()})
                break

    # Re-apply the expected capitalization for downstream lookups
    df = df.rename(columns={
        "a54 count": "A54 count",
        "u54 count": "U54 count",
        "total ala agc count": "total ala-agc count",
    })

    # Keep only rows with known lineage labels and valid taxids
    df = df[df["lineage"].notna()].copy()
    df["lineage"] = df["lineage"].astype(str).str.strip()
    df = df[df["lineage"].str.upper() != "N/A"]
    print(f"  {len(df)} species after removing N/A entries")

    df["ncbi taxid"] = pd.to_numeric(df["ncbi taxid"], errors="coerce")
    df = df.dropna(subset=["ncbi taxid"]).copy()
    df["ncbi taxid"] = df["ncbi taxid"].astype(int)
    df = df.drop_duplicates(subset=["ncbi taxid"]).copy()
 
    # --- NCBI taxonomy DB ------------------------------------------------------
    from ete3 import NCBITaxa
    ncbi = NCBITaxa()
 
    if args.update_db:
        print("Updating NCBI taxonomy database (this may take a few minutes)...")
        ncbi.update_taxonomy_database()
 
    # --- Build tree ------------------------------------------------------
    taxids = df["ncbi taxid"].tolist()
    tree = build_tree_from_taxids(taxids)
 
    # --- Style & render ------------------------------------------------------
    ts = style_tree(tree, df)
    ts.mode = args.mode
 
    out = args.output

    print(f"Rendering tree to {out}...")
    tree.render(out, tree_style=ts, dpi=600, units="px", w=2400)
    print(f"Done! Tree saved to: {out}")
 
 
if __name__ == "__main__":
    main()