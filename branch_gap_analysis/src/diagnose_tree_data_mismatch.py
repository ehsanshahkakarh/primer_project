#!/usr/bin/env python3
"""
Diagnose the mismatch between tree colorstrips and actual data.

This script compares:
1. What taxa are in the tree (genera)
2. What taxa are in the final_merger files (aggregated at different levels)
3. The matching logic in create_itol_annotations.py

Goal: Identify why division/family NF colors don't match actual genus-level data.
"""

import pandas as pd
import csv
from pathlib import Path
from collections import defaultdict

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
GENUS_PARSE_FILE = REPO_ROOT / "00_gaps_taxonomic/00parse_database/eukcensus_parse/18S_censusparse/output/eukcensus_18S_by_genus.csv"
FINAL_MERGER_BASE = REPO_ROOT / "00_gaps_taxonomic/00parse_database/final_merger/outputs"

DIVISION_RANKS = {'phylum', 'division', 'clade', 'kingdom', 'subkingdom'}
FAMILY_RANKS = {'family', 'subfamily', 'superfamily', 'tribe', 'subtribe'}

def main():
    print("=" * 80)
    print("DIAGNOSIS: Tree Colorstrip vs Final Merger Data Mismatch")
    print("=" * 80)
    
    # Load final_merger data
    div_df = pd.read_csv(FINAL_MERGER_BASE / "18s_ncbi_merged_division.csv")
    fam_df = pd.read_csv(FINAL_MERGER_BASE / "18s_ncbi_merged_family.csv")
    gen_df = pd.read_csv(FINAL_MERGER_BASE / "18s_ncbi_merged_genus.csv")
    
    print(f"\nFinal Merger Data Summary:")
    print(f"  Division level: {len(div_df)} unique taxa")
    print(f"  Family level: {len(fam_df)} unique taxa")
    print(f"  Genus level: {len(gen_df)} unique taxa")
    
    # Create lookup sets
    div_lookup = set(div_df['division'].values)
    fam_lookup = set(fam_df['family'].values)
    gen_lookup = set(gen_df['genus'].values)
    
    # Load census data to analyze lineages
    census_df = pd.read_csv(GENUS_PARSE_FILE)
    print(f"\nCensus Data:")
    print(f"  Total genera in tree: {len(census_df)}")
    
    # Analyze matching for sample genera
    sample_genera = ['Acanthamoeba', 'Amphitrema', 'Arcellinida.U.genus', 'Amoebozoa.U.genus']
    
    print(f"\n{'='*80}")
    print("SAMPLE GENUS ANALYSIS")
    print(f"{'='*80}")
    
    for genus_name in sample_genera:
        genus_rows = census_df[census_df['Name_to_use'] == genus_name]
        if genus_rows.empty:
            print(f"\n{genus_name}: NOT FOUND IN CENSUS")
            continue
        
        row = genus_rows.iloc[0]
        lineage = row['lineage']
        lineage_ranks = row['lineage_ranks']
        
        names = [x.strip() for x in lineage.split(';')]
        ranks = [x.strip() for x in lineage_ranks.split(';')]
        
        print(f"\n{genus_name}:")
        print(f"  Lineage: {lineage}")
        print(f"  Ranks: {lineage_ranks}")
        
        # Check division matching
        divisions_in_lineage = []
        for name, rank in zip(names, ranks):
            if rank.lower() in DIVISION_RANKS:
                divisions_in_lineage.append(name)
        
        div_match = None
        for div in divisions_in_lineage:
            if div in div_lookup:
                div_match = div
                break
        
        print(f"  Division-level taxa in lineage: {divisions_in_lineage}")
        print(f"  Division match in final_merger: {div_match}")
        if div_match:
            div_nf = div_df[div_df['division'] == div_match]['novelty_factor'].iloc[0]
            print(f"    NF = {div_nf}")
        
        # Check family matching
        families_in_lineage = []
        for name, rank in zip(names, ranks):
            if rank.lower() in FAMILY_RANKS:
                families_in_lineage.append(name)
        
        fam_match = None
        for fam in families_in_lineage:
            if fam in fam_lookup:
                fam_match = fam
                break
        
        print(f"  Family-level taxa in lineage: {families_in_lineage}")
        print(f"  Family match in final_merger: {fam_match}")
        if fam_match:
            fam_nf = fam_df[fam_df['family'] == fam_match]['novelty_factor'].iloc[0]
            print(f"    NF = {fam_nf}")
        
        # Check genus matching
        genus_match = genus_name in gen_lookup
        print(f"  Genus match in final_merger: {genus_match}")
        if genus_match:
            gen_nf = gen_df[gen_df['genus'] == genus_name]['novelty_factor'].iloc[0]
            print(f"    NF = {gen_nf}")
    
    # Overall statistics
    print(f"\n{'='*80}")
    print("OVERALL MATCHING STATISTICS")
    print(f"{'='*80}")
    
    div_match_count = 0
    fam_match_count = 0
    gen_match_count = 0
    total_genera = 0
    
    for idx, row in census_df.iterrows():
        total_genera += 1
        genus_name = row['Name_to_use']
        lineage = row['lineage']
        lineage_ranks = row['lineage_ranks']
        
        names = [x.strip() for x in lineage.split(';')]
        ranks = [x.strip() for x in lineage_ranks.split(';')]
        
        # Check division
        for name, rank in zip(names, ranks):
            if rank.lower() in DIVISION_RANKS and name in div_lookup:
                div_match_count += 1
                break
        
        # Check family
        for name, rank in zip(names, ranks):
            if rank.lower() in FAMILY_RANKS and name in fam_lookup:
                fam_match_count += 1
                break
        
        # Check genus
        if genus_name in gen_lookup:
            gen_match_count += 1
    
    print(f"Total genera in census: {total_genera}")
    print(f"Genera with division-level NF data: {div_match_count} ({div_match_count/total_genera*100:.1f}%)")
    print(f"Genera with family-level NF data: {fam_match_count} ({fam_match_count/total_genera*100:.1f}%)")
    print(f"Genera with genus-level NF data: {gen_match_count} ({gen_match_count/total_genera*100:.1f}%)")
    
    print(f"\n{'='*80}")
    print("CONCLUSION:")
    print(f"{'='*80}")
    print(f"The final_merger files aggregate data at HIGHER taxonomic levels,")
    print(f"while the tree contains MANY MORE genera. This explains the disconnect:")
    print(f"  - Division colorstrip: Only {div_match_count}/{total_genera} genera get colored")
    print(f"  - Family colorstrip: Only {fam_match_count}/{total_genera} genera get colored")
    print(f"  - Genus colorstrip: Only {gen_match_count}/{total_genera} genera get colored")

if __name__ == "__main__":
    main()
