#!/usr/bin/env python3
"""
Create filtered annotation files for each domain (Archaea, Bacteria).
"""

import re
from datetime import datetime
from pathlib import Path

def get_taxa_from_tree(tree_file):
    """Extract all taxon names from a newick tree file."""
    with open(tree_file, 'r') as f:
        tree_text = f.read()
    # Match taxon names with rank suffixes like Name_G, Name_F, etc.
    pattern = r'([A-Za-z0-9_/]+_[A-Z])[,);]'
    taxa = set(re.findall(pattern, tree_text))
    return taxa

def filter_annotation_file(input_file, output_file, valid_taxa, domain_name):
    """Filter an annotation file to only include taxa in the tree."""
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Identify header lines (before DATA)
    header_lines = []
    data_started = False
    data_lines = []
    
    for line in lines:
        if line.strip() == 'DATA':
            header_lines.append(line)
            data_started = True
        elif not data_started:
            header_lines.append(line)
        else:
            data_lines.append(line)
    
    # Filter data lines - keep only those with taxa in our tree
    filtered_data = []
    for line in data_lines:
        if not line.strip():
            continue
        parts = line.split('\t')
        if parts:
            taxon = parts[0].strip()
            if taxon in valid_taxa:
                filtered_data.append(line)
    
    # Write output
    with open(output_file, 'w') as f:
        for line in header_lines:
            f.write(line)
        for line in filtered_data:
            f.write(line)
    
    return len(filtered_data)

def main():
    script_dir = Path(__file__).parent
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Annotation files to process
    annotation_files = [
        'full_taxonomy_colors.txt',
        'full_taxonomy_colorstrip.txt', 
        'full_taxonomy_internal_labels.txt',
        'full_taxonomy_labels.txt',
    ]
    
    domains = ['archaea', 'bacteria']
    
    for domain in domains:
        domain_dir = script_dir / domain
        
        # Find the tree file
        tree_files = list(domain_dir.glob(f"{domain}_tree_*.nwk"))
        if not tree_files:
            print(f"No tree file found for {domain}")
            continue
        
        tree_file = tree_files[0]
        print(f"\n=== Processing {domain.upper()} ===")
        print(f"Tree file: {tree_file.name}")
        
        # Get taxa from tree
        taxa = get_taxa_from_tree(tree_file)
        print(f"Taxa in tree: {len(taxa)}")
        
        # Filter each annotation file
        for ann_file in annotation_files:
            input_path = script_dir / ann_file
            if not input_path.exists():
                print(f"  Skipping {ann_file} (not found)")
                continue
            
            # Create output filename
            output_name = ann_file.replace('full_taxonomy_', f'{domain}_')
            output_path = domain_dir / output_name
            
            count = filter_annotation_file(input_path, output_path, taxa, domain)
            print(f"  {output_name}: {count} entries")

if __name__ == "__main__":
    main()

