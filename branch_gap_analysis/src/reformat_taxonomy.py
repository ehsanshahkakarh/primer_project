#!/usr/bin/env python3
"""
Reformat taxonomy data with separate columns for each taxonomic rank.

Uses lineage data from the census files (already parsed) and joins with
novelty factors from the merged files. Outputs CSVs with individual columns
for each taxonomic rank found in the data.

Output: output/{16s,18s}/taxa_with_ranks.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
from config import (
    SOURCE_16S_FILES, SOURCE_18S_FILES,
    METADATA_DIR, OUTPUT_DIR, OTU_PIPELINE_ROOT
)

# Census file paths (note: 16S uses 'eukcensus16S_', 18S uses 'eukcensus_18S_')
CENSUS_16S_FILES = {
    "division": OTU_PIPELINE_ROOT / "eukcensus_parse" / "16S_censusparse" / "output" / "eukcensus16S_by_division.csv",
    "family": OTU_PIPELINE_ROOT / "eukcensus_parse" / "16S_censusparse" / "output" / "eukcensus16S_by_family.csv",
    "genus": OTU_PIPELINE_ROOT / "eukcensus_parse" / "16S_censusparse" / "output" / "eukcensus16S_by_genus.csv",
}

CENSUS_18S_FILES = {
    "division": OTU_PIPELINE_ROOT / "eukcensus_parse" / "18S_censusparse" / "output" / "eukcensus_18S_by_division.csv",
    "family": OTU_PIPELINE_ROOT / "eukcensus_parse" / "18S_censusparse" / "output" / "eukcensus_18S_by_family.csv",
    "genus": OTU_PIPELINE_ROOT / "eukcensus_parse" / "18S_censusparse" / "output" / "eukcensus_18S_by_genus.csv",
}


def parse_census_lineage(lineage_str: str, lineage_ranks_str: str) -> dict[str, str]:
    """Parse census lineage and ranks into a dictionary.

    Census format:
        lineage: "cellular organisms;Eukaryota;Opisthokonta;Metazoa"
        lineage_ranks: "cellular root;domain;clade;kingdom"

    Returns dict mapping rank -> taxon name for all ranks found.
    """
    if pd.isna(lineage_str) or pd.isna(lineage_ranks_str):
        return {}

    names = lineage_str.split(";")
    ranks = lineage_ranks_str.split(";")

    rank_dict = {}
    for name, rank in zip(names, ranks):
        name = name.strip()
        rank = rank.strip().lower()
        if name and rank and rank != "original_name":
            # Normalize some rank names
            if rank == "cellular root":
                continue  # Skip this
            if rank == "superkingdom":
                rank = "domain"
            rank_dict[rank] = name

    return rank_dict


def load_census_lineages(census_files: dict) -> dict[str, dict[str, str]]:
    """Load lineage data from census files.

    Returns dict mapping taxon_name -> rank_dict
    """
    name_to_lineage = {}

    for level, filepath in census_files.items():
        if filepath.exists():
            df = pd.read_csv(filepath)
            if 'Name_to_use' in df.columns and 'lineage' in df.columns and 'lineage_ranks' in df.columns:
                for _, row in df.iterrows():
                    name = row['Name_to_use']
                    lineage = row['lineage']
                    ranks = row['lineage_ranks']
                    rank_dict = parse_census_lineage(lineage, ranks)
                    if rank_dict:
                        name_to_lineage[name] = rank_dict
            print(f"    Loaded lineages from {level}: {len(df)} entries")

    return name_to_lineage


def is_unresolved_taxon(name: str) -> bool:
    """Check if taxon name is unresolved/placeholder (should be filtered out).

    Filters:
    - .U. pattern (e.g., 'Euglenida.U.family')
    - _X or _XX suffix patterns (e.g., 'Monothalamids_X', 'Something_XX')
    """
    import re
    if pd.isna(name):
        return True
    # Filter .U. pattern
    if '.U.' in name:
        return True
    # Filter _X, _XX, _XXX etc. suffix patterns
    if re.search(r'_X+$', name):
        return True
    return False


def load_matched_taxa(source_files: dict, filter_unresolved: bool = True) -> pd.DataFrame:
    """Load matched taxa from merged files with novelty factors and counts.

    Args:
        source_files: Dictionary mapping level -> filepath
        filter_unresolved: If True, filter out .U. and _XX taxa
    """
    all_taxa = []

    # Columns to extract (with fallbacks if not present)
    count_columns = [
        'census_otu_count',
        'census_size_count',
        'ncbi_genome_count',
        'ncbi_species_count'
    ]

    for level, filepath in source_files.items():
        if filepath.exists():
            df = pd.read_csv(filepath)
            taxon_col = level
            if taxon_col in df.columns:
                # Filter for matched OR census_only taxa (include inf NF)
                if 'match_status' in df.columns:
                    df = df[df['match_status'].isin(['matched', 'census_only'])]

                # Build list of columns to extract
                cols_to_extract = [taxon_col, 'novelty_factor', 'domain', 'census_taxid']
                for col in count_columns:
                    if col in df.columns:
                        cols_to_extract.append(col)

                subset = df[cols_to_extract].copy()
                subset = subset.rename(columns={taxon_col: 'taxon_name'})
                subset['level'] = level
                # Keep infinite novelty factors (census-only taxa) but filter out NaN
                # Replace inf with a large number for sorting/visualization
                subset = subset[subset['novelty_factor'].notna()]
                # Convert inf to a high value (9999) for processing
                subset['novelty_factor'] = subset['novelty_factor'].replace([np.inf], 9999.0)

                # Filter out unresolved taxa
                if filter_unresolved:
                    before_count = len(subset)
                    subset = subset[~subset['taxon_name'].apply(is_unresolved_taxon)]
                    filtered_count = before_count - len(subset)
                    if filtered_count > 0:
                        print(f"    Filtered {filtered_count} unresolved taxa from {level}")

                all_taxa.append(subset)
                print(f"    Loaded {len(subset)} matched taxa from {level}")

    if all_taxa:
        return pd.concat(all_taxa, ignore_index=True)
    return pd.DataFrame()


def collect_all_ranks(name_to_lineage: dict[str, dict[str, str]]) -> list[str]:
    """Collect all unique ranks found across all lineages, in a sensible order."""
    # Preferred rank order
    rank_order = [
        "domain", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum",
        "infraphylum", "superclass", "class", "subclass", "infraclass",
        "superorder", "order", "suborder", "infraorder", "parvorder",
        "superfamily", "family", "subfamily", "tribe", "subtribe",
        "genus", "subgenus", "species_group", "species_subgroup", "species",
        "subspecies", "strain", "clade", "no rank"
    ]

    # Collect all ranks found
    all_ranks = set()
    for rank_dict in name_to_lineage.values():
        all_ranks.update(rank_dict.keys())

    # Sort by preferred order, then alphabetically for unknown ranks
    def rank_sort_key(r):
        if r in rank_order:
            return (0, rank_order.index(r))
        return (1, r)

    return sorted(all_ranks, key=rank_sort_key)


def process_dataset(name: str, source_files: dict, census_files: dict, output_dir: Path):
    """Process a dataset and create reformatted taxonomy CSV."""
    print(f"\n{'='*70}")
    print(f"Processing {name}")
    print(f"{'='*70}")

    # Load matched taxa with novelty factors
    taxa_df = load_matched_taxa(source_files)
    if taxa_df.empty:
        print("  No taxa loaded!")
        return
    print(f"  Total: {len(taxa_df)} matched taxa loaded")

    # Load lineages from census files
    print("  Loading census lineages...")
    name_to_lineage = load_census_lineages(census_files)
    print(f"  Loaded {len(name_to_lineage)} lineages from census")

    # Collect all ranks found
    all_ranks = collect_all_ranks(name_to_lineage)
    print(f"  Found {len(all_ranks)} unique ranks: {all_ranks[:10]}...")

    # Create output dataframe with rank columns
    output_rows = []
    matched_lineages = 0

    for _, row in taxa_df.iterrows():
        taxon_name = row['taxon_name']
        rank_dict = name_to_lineage.get(taxon_name, {})

        if rank_dict:
            matched_lineages += 1

        # Convert taxid to int if it's a valid number
        taxid = row['census_taxid']
        if pd.notna(taxid):
            taxid = int(taxid)

        out_row = {
            'taxon_name': taxon_name,
            'level': row['level'],
            'novelty_factor': row['novelty_factor'],
            'source_domain': row['domain'],
            'census_taxid': taxid,
            # Census counts
            'census_otu_count': row.get('census_otu_count', pd.NA),
            'census_size_count': row.get('census_size_count', pd.NA),
            # NCBI counts
            'ncbi_genome_count': row.get('ncbi_genome_count', pd.NA),
            'ncbi_species_count': row.get('ncbi_species_count', pd.NA),
        }

        # Add each rank as a column
        for rank in all_ranks:
            out_row[rank] = rank_dict.get(rank, pd.NA)

        output_rows.append(out_row)

    print(f"  Matched lineages for {matched_lineages}/{len(taxa_df)} taxa")

    # Create output DataFrame
    output_df = pd.DataFrame(output_rows)

    # Sort by species_count descending, then novelty_factor descending
    # This highlights taxa with many NCBI species AND high underrepresentation
    output_df = output_df.sort_values(
        ['ncbi_species_count', 'novelty_factor'],
        ascending=[False, False]
    )

    # Save to output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save as CSV
    csv_file = output_dir / "taxa_with_ranks.csv"
    output_df.to_csv(csv_file, index=False)
    print(f"  Saved CSV: {csv_file}")

    # Save as TSV
    tsv_file = output_dir / "taxa_with_ranks.tsv"
    output_df.to_csv(tsv_file, index=False, sep='\t')
    print(f"  Saved TSV: {tsv_file}")

    print(f"  Columns: {list(output_df.columns)}")

    # Also save a summary by domain
    if 'source_domain' in output_df.columns:
        domain_summary = output_df.groupby('source_domain').agg({
            'taxon_name': 'count',
            'novelty_factor': ['mean', 'max']
        }).round(2)
        summary_file = output_dir / "domain_summary.csv"
        domain_summary.to_csv(summary_file)
        print(f"  Saved domain summary: {summary_file}")


def main():
    # Process 16S
    process_dataset(
        "16S",
        SOURCE_16S_FILES,
        CENSUS_16S_FILES,
        METADATA_DIR / "16s"
    )

    # Process 18S
    process_dataset(
        "18S",
        SOURCE_18S_FILES,
        CENSUS_18S_FILES,
        METADATA_DIR / "18s"
    )

    print(f"\n{'='*70}")
    print("Taxonomy reformatting complete!")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()

