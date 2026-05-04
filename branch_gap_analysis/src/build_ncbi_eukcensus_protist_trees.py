#!/usr/bin/env python3
"""
Build two genus-level Newick trees for 18S work:

1) NCBI protists — from ncbi_parse species_grouped CSV. A row counts as protist if
   lineage_taxids include Eukaryota (2759) and exclude Metazoa (33208) and Fungi (4751).
   Each species row contributes its lineage truncated at rank "genus".

2) EuKCensus — from eukcensus_18S_by_genus.csv (same logic as build_genus_tree.py).

Outputs (under branch_gap_analysis/output/18s/tree/):
  - ncbi_protist_genus_tree.nwk
  - eukcensus_18s_genus_tree.nwk
"""

from __future__ import annotations

import csv
from pathlib import Path

from build_genus_tree import (
    TaxonNode,
    count_nodes,
    integrate_genus_lineage,
    print_tree_summary,
    should_skip_taxon,
    tree_to_newick,
)


# NCBI taxids: Eukaryota, Metazoa, Fungi (operational "protist" = Eukaryota \\ (Metazoa ∪ Fungi))
NCBI_EUKARYOTA = "2759"
NCBI_METAZOA = "33208"
NCBI_FUNGI = "4751"

DEFAULT_SPECIES_GROUPED = (
    "species_grouped_20260301_214126.csv"
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent.parent.parent


def ncbi_lineage_is_protist(lineage_taxids: str) -> bool:
    ids = {x.strip() for x in (lineage_taxids or "").split(";") if x.strip()}
    if NCBI_EUKARYOTA not in ids:
        return False
    if NCBI_METAZOA in ids or NCBI_FUNGI in ids:
        return False
    return True


def truncate_to_genus(
    lineage: str,
    lineage_ranks: str,
    lineage_taxids: str,
) -> tuple[list[str], list[str], list[str]] | None:
    names = lineage.split(";")
    ranks = lineage_ranks.split(";")
    taxids = (
        lineage_taxids.split(";")
        if lineage_taxids
        else [""] * len(names)
    )
    if len(names) != len(ranks):
        return None
    if len(taxids) < len(names):
        taxids = taxids + [""] * (len(names) - len(taxids))
    genus_idx = None
    for i, r in enumerate(ranks):
        if r.strip() == "genus":
            genus_idx = i
    if genus_idx is None:
        return None
    return (
        names[: genus_idx + 1],
        ranks[: genus_idx + 1],
        taxids[: genus_idx + 1],
    )


def build_ncbi_protist_tree(species_csv: Path) -> tuple[TaxonNode, dict]:
    root = TaxonNode(name="Life", rank="root")
    stats = {"total": 0, "skipped": 0, "included": 0}

    with species_csv.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            stats["total"] += 1
            lt = row.get("lineage_taxids", "") or ""
            if not ncbi_lineage_is_protist(lt):
                stats["skipped"] += 1
                continue

            lineage = row.get("lineage", "") or ""
            lr = row.get("lineage_ranks", "") or ""
            truncated = truncate_to_genus(lineage, lr, lt)
            if truncated is None:
                stats["skipped"] += 1
                continue

            names, ranks, taxids = truncated
            genus_name = names[-1].strip() if names else ""
            if not genus_name or should_skip_taxon(genus_name):
                stats["skipped"] += 1
                continue

            otu = int(row.get("total_genome_count", 0) or 0)
            sz = int(row.get("isolate_genome_count", 0) or 0)
            integrate_genus_lineage(root, names, ranks, taxids, otu, sz)
            stats["included"] += 1

    return root, stats


def build_eukcensus_tree(eukcensus_csv: Path) -> tuple[TaxonNode, dict]:
    root = TaxonNode(name="Life", rank="root")
    stats = {"total": 0, "skipped": 0, "included": 0}

    with eukcensus_csv.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            stats["total"] += 1
            genus_name = row.get("Name_to_use", "")
            lineage = row.get("lineage", "")
            lineage_ranks = row.get("lineage_ranks", "")
            lineage_taxids = row.get("lineage_taxids", "")
            otu_count = int(row.get("otu_count", 0) or 0)
            size_count = int(row.get("size_count", 0) or 0)

            if not lineage or not lineage_ranks:
                stats["skipped"] += 1
                continue
            if should_skip_taxon(genus_name):
                stats["skipped"] += 1
                continue

            stats["included"] += 1
            names = lineage.split(";")
            ranks = lineage_ranks.split(";")
            taxids = (
                lineage_taxids.split(";")
                if lineage_taxids
                else [""] * len(names)
            )
            integrate_genus_lineage(root, names, ranks, taxids, otu_count, size_count)

    return root, stats


def main() -> None:
    root_data = _repo_root()
    ncbi_csv = (
        root_data
        / "00_gaps_taxonomic/00parse_database/ncbi_parse/output"
        / DEFAULT_SPECIES_GROUPED
    )
    euk_csv = (
        root_data
        / "00_gaps_taxonomic/00parse_database/eukcensus_parse/18S_censusparse/output"
        / "eukcensus_18S_by_genus.csv"
    )
    out_dir = (
        root_data
        / "primer_project/branch_gap_analysis/output/18s/tree"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("Tree 1: NCBI protist genera (!= Metazoa, != Fungi; Eukaryota only)")
    print("=" * 72)
    print(f"Reading: {ncbi_csv}")
    ncbi_root, nstats = build_ncbi_protist_tree(ncbi_csv)
    print(
        f"  Rows: {nstats['total']} | skipped: {nstats['skipped']} | integrated: {nstats['included']}"
    )
    nt, nterm = count_nodes(ncbi_root)
    print(f"  Nodes: {nt} total, {nterm} terminal (genus) tips")
    print_tree_summary(ncbi_root, max_depth=2)
    ncbi_nwk = out_dir / "ncbi_protist_genus_tree.nwk"
    ncbi_nwk.write_text(tree_to_newick(ncbi_root) + ";", encoding="utf-8")
    print(f"Saved: {ncbi_nwk}\n")

    print("=" * 72)
    print("Tree 2: EuKCensus 18S genera")
    print("=" * 72)
    print(f"Reading: {euk_csv}")
    euk_root, estats = build_eukcensus_tree(euk_csv)
    print(
        f"  Rows: {estats['total']} | skipped: {estats['skipped']} | integrated: {estats['included']}"
    )
    et, eterm = count_nodes(euk_root)
    print(f"  Nodes: {et} total, {eterm} terminal (genus) tips")
    print_tree_summary(euk_root, max_depth=2)
    euk_nwk = out_dir / "eukcensus_18s_genus_tree.nwk"
    euk_nwk.write_text(tree_to_newick(euk_root) + ";", encoding="utf-8")
    print(f"Saved: {euk_nwk}")


if __name__ == "__main__":
    main()
