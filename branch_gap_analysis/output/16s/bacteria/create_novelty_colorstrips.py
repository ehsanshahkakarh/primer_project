#!/usr/bin/env python3
"""
Create separate novelty colorstrip files for Phylum, Family, and Genus ranks.
Uses the orange-red-blue color scheme based on novelty values.
"""

import re
from pathlib import Path

def get_color_for_novelty(nf_value):
    """
    Return color based on novelty factor value.
    Uses the same scheme as build_full_taxonomy_tree.py:
    - Gray (#888888): NF <= 1 (normal representation)
    - Orange (#FF8000): NF 1-2 (up to 100% underrepresented)
    - Orange-Red gradient: NF 2-11 (100% to 1000% underrepresented)
    - Navy Blue (#000080): NF > 11 (egregiously underrepresented)
    """
    if nf_value <= 1:
        return "#888888"  # Gray - normal representation
    elif nf_value <= 2:
        return "#FF8000"  # Orange
    elif nf_value <= 11:
        # Red gradient: NF 2 to 11 (100% to 1000% underrepresented)
        # Gradient from orange-red to dark red
        norm = (nf_value - 2) / (11 - 2)  # 0 to 1
        r = int(255 - (55 * norm))  # 255 -> 200
        g = int(80 * (1 - norm))    # 80 -> 0
        return f"#{r:02x}{g:02x}00"
    else:
        return "#000080"  # Navy blue - extreme novelty (>1000%)

def get_taxa_from_tree(tree_file):
    """Extract all taxon names from a newick tree file."""
    with open(tree_file, 'r') as f:
        tree_text = f.read()
    pattern = r'([A-Za-z0-9_/]+_[A-Z])[,);]'
    return set(re.findall(pattern, tree_text))

def parse_labels_file(labels_file):
    """Parse the labels file to extract novelty values."""
    novelty_data = {}
    with open(labels_file, 'r') as f:
        for line in f:
            if 'NF=' in line:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    taxon = parts[0]
                    # Extract NF value
                    match = re.search(r'NF=([0-9.]+)', parts[1])
                    if match:
                        novelty_data[taxon] = float(match.group(1))
    return novelty_data

def create_colorstrip(output_file, taxa_data, label, strip_width=25):
    """Create a colorstrip file for iTOL."""
    with open(output_file, 'w') as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR TAB\n")
        f.write(f"DATASET_LABEL\t{label}\n")
        f.write("COLOR\t#cc0000\n")
        f.write(f"STRIP_WIDTH\t{strip_width}\n")
        f.write("MARGIN\t5\n")
        f.write("SHOW_INTERNAL\t1\n")  # Show on internal nodes too
        f.write("DATA\n")
        
        for taxon, nf_value in sorted(taxa_data.items()):
            color = get_color_for_novelty(nf_value)
            f.write(f"{taxon}\t{color}\t{nf_value}\n")
    
    return len(taxa_data)

def main():
    script_dir = Path(__file__).parent
    parent_dir = script_dir.parent
    
    # Find tree file
    tree_files = list(script_dir.glob("bacteria_tree_*.nwk"))
    if not tree_files:
        print("No bacteria tree file found!")
        return
    
    tree_file = tree_files[0]
    print(f"Using tree: {tree_file.name}")
    
    # Get taxa in tree
    tree_taxa = get_taxa_from_tree(tree_file)
    print(f"Taxa in tree: {len(tree_taxa)}")
    
    # Parse labels file for novelty values
    labels_file = parent_dir / "full_taxonomy_labels.txt"
    all_novelty = parse_labels_file(labels_file)
    print(f"Total novelty entries: {len(all_novelty)}")
    
    # Filter by rank and tree membership
    phylum_data = {k: v for k, v in all_novelty.items() if k.endswith('_P') and k in tree_taxa}
    family_data = {k: v for k, v in all_novelty.items() if k.endswith('_F') and k in tree_taxa}
    genus_data = {k: v for k, v in all_novelty.items() if k.endswith('_G') and k in tree_taxa}
    
    print(f"\nFiltered for bacteria tree:")
    print(f"  Phyla: {len(phylum_data)}")
    print(f"  Families: {len(family_data)}")
    print(f"  Genera: {len(genus_data)}")
    
    # Create colorstrip files
    print("\nCreating colorstrip files...")
    
    count = create_colorstrip(
        script_dir / "phylum_novelty_colorstrip.txt",
        phylum_data,
        "Phylum Novelty",
        strip_width=30
    )
    print(f"  phylum_novelty_colorstrip.txt: {count} entries")
    
    count = create_colorstrip(
        script_dir / "family_novelty_colorstrip.txt",
        family_data,
        "Family Novelty",
        strip_width=25
    )
    print(f"  family_novelty_colorstrip.txt: {count} entries")
    
    count = create_colorstrip(
        script_dir / "genus_novelty_colorstrip.txt",
        genus_data,
        "Genus Novelty",
        strip_width=20
    )
    print(f"  genus_novelty_colorstrip.txt: {count} entries")
    
    print("\nDone! Load files in iTOL in this order for inside-to-outside rings:")
    print("  1. phylum_novelty_colorstrip.txt (innermost)")
    print("  2. family_novelty_colorstrip.txt (middle)")
    print("  3. genus_novelty_colorstrip.txt (outermost)")

if __name__ == "__main__":
    main()

