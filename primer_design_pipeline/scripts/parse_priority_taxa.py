#!/usr/bin/env python3
"""
Parse Priority Taxa from iTOL Prune Lists

This script parses the dark_blue priority taxa from iTOL prune list files and
extracts taxonomic information based on rank suffixes:
  _P = Phylum/Division (maps to 'division' column in cluster TSV)
  _C = Class (not directly in cluster TSV - need to query by lineage)
  _O = Order (not directly in cluster TSV - need to query by lineage)
  _F = Family (maps to 'family' column)
  _G = Genus (maps to 'genus' column)
  _S = Species (subset of genus)

Usage:
    python parse_priority_taxa.py <prune_list_file> <output_json>
"""

import sys
import json
import re
from pathlib import Path

# Rank suffix to column/query type mapping
RANK_SUFFIX_MAP = {
    '_P': {'rank': 'phylum', 'query_column': 'division', 'query_type': 'direct'},
    '_C': {'rank': 'class', 'query_column': None, 'query_type': 'lineage'},
    '_O': {'rank': 'order', 'query_column': None, 'query_type': 'lineage'},
    '_F': {'rank': 'family', 'query_column': 'family', 'query_type': 'direct'},
    '_G': {'rank': 'genus', 'query_column': 'genus', 'query_type': 'direct'},
    '_S': {'rank': 'species', 'query_column': 'genus', 'query_type': 'partial'},
    '_I': {'rank': 'infraclass', 'query_column': None, 'query_type': 'lineage'},  # e.g., Teleostei_I
}

def parse_rank_suffix(taxon_name: str) -> dict:
    """Parse a taxon name to extract base name and rank information."""
    for suffix, info in RANK_SUFFIX_MAP.items():
        if taxon_name.endswith(suffix):
            base_name = taxon_name[:-len(suffix)]
            return {
                'original_name': taxon_name,
                'base_name': base_name,
                'suffix': suffix,
                **info
            }
    
    # No recognized suffix - treat as unknown
    return {
        'original_name': taxon_name,
        'base_name': taxon_name,
        'suffix': None,
        'rank': 'unknown',
        'query_column': None,
        'query_type': 'unknown'
    }

def parse_prune_list(prune_list_file: str) -> list:
    """Parse an iTOL prune list file and extract taxa with rank info."""
    taxa = []
    
    with open(prune_list_file, 'r') as f:
        in_data_section = False
        for line in f:
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
            
            # Look for DATA section marker
            if line == 'DATA':
                in_data_section = True
                continue
            
            # Skip header lines before DATA
            if not in_data_section:
                continue
            
            # Parse taxon name
            taxon_info = parse_rank_suffix(line)
            taxa.append(taxon_info)
    
    return taxa

def summarize_taxa(taxa: list) -> dict:
    """Generate summary statistics for parsed taxa."""
    summary = {'total': len(taxa), 'by_rank': {}, 'by_query_type': {}}
    
    for taxon in taxa:
        rank = taxon['rank']
        qtype = taxon['query_type']
        
        summary['by_rank'][rank] = summary['by_rank'].get(rank, 0) + 1
        summary['by_query_type'][qtype] = summary['by_query_type'].get(qtype, 0) + 1
    
    return summary

def main():
    if len(sys.argv) < 2:
        print("Usage: python parse_priority_taxa.py <prune_list_file> [output_json]")
        sys.exit(1)
    
    prune_list_file = sys.argv[1]
    output_json = sys.argv[2] if len(sys.argv) > 2 else None
    
    print(f"📖 Parsing priority taxa from: {prune_list_file}")
    
    taxa = parse_prune_list(prune_list_file)
    summary = summarize_taxa(taxa)
    
    print(f"\n📊 Summary:")
    print(f"   Total taxa: {summary['total']}")
    print(f"\n   By rank:")
    for rank, count in sorted(summary['by_rank'].items()):
        print(f"     {rank}: {count}")
    print(f"\n   By query type:")
    for qtype, count in sorted(summary['by_query_type'].items()):
        print(f"     {qtype}: {count}")
    
    # Output JSON if requested
    if output_json:
        output = {'taxa': taxa, 'summary': summary}
        with open(output_json, 'w') as f:
            json.dump(output, f, indent=2)
        print(f"\n💾 Saved to: {output_json}")
    
    # Also print first few taxa as preview
    print(f"\n📋 First 10 taxa:")
    for taxon in taxa[:10]:
        print(f"   {taxon['original_name']} → {taxon['base_name']} ({taxon['rank']}, {taxon['query_type']})")
    
    return taxa

if __name__ == "__main__":
    main()

