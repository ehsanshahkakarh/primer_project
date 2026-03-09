#!/usr/bin/env python3
"""
Split the 16S full taxonomy tree into domain-specific subtrees (Archaea, Bacteria).
"""

import re
from datetime import datetime
from pathlib import Path

def extract_domain_subtree(newick_text, domain_name):
    """Extract the subtree for a specific domain from the newick tree."""
    # Pattern matches )Domain_D at the end of a clade
    pattern = rf'\){re.escape(domain_name)}_D'
    match = re.search(pattern, newick_text)
    
    if not match:
        return None
    
    end_pos = match.start()
    paren_count = 1
    start_pos = end_pos - 1
    
    while start_pos >= 0 and paren_count > 0:
        if newick_text[start_pos] == ')':
            paren_count += 1
        elif newick_text[start_pos] == '(':
            paren_count -= 1
        start_pos -= 1
    
    start_pos += 1
    subtree_content = newick_text[start_pos:end_pos + 1]
    subtree = subtree_content + f"{domain_name}_D;"
    
    return subtree

def count_leaves(newick_text, suffix="_G"):
    """Count number of leaves (genera by default) in a newick tree."""
    pattern = rf'[A-Za-z0-9_]+{re.escape(suffix)}[,)]'
    return len(re.findall(pattern, newick_text))

def main():
    script_dir = Path(__file__).parent
    
    # Find the tree file
    tree_file = script_dir / "full_taxonomy_tree.nwk"
    if not tree_file.exists():
        print("Error: full_taxonomy_tree.nwk not found!")
        return
    
    print(f"Processing: {tree_file.name}")
    
    with open(tree_file, 'r') as f:
        newick_text = f.read().strip()
    
    print(f"Original tree size: {len(newick_text):,} characters")
    print(f"Original genera count: {count_leaves(newick_text)}")
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Extract Archaea
    print("\nExtracting Archaea...")
    archaea_subtree = extract_domain_subtree(newick_text, "Archaea")
    if archaea_subtree:
        archaea_dir = script_dir / "archaea"
        archaea_dir.mkdir(exist_ok=True)
        output_file = archaea_dir / f"archaea_tree_{timestamp}.nwk"
        with open(output_file, 'w') as f:
            f.write(archaea_subtree)
        print(f"  Created: {output_file.name}")
        print(f"  Genera: {count_leaves(archaea_subtree)}")
        print(f"  Size: {len(archaea_subtree):,} characters")
    else:
        print("  ERROR: Could not extract Archaea subtree")
    
    # Extract Bacteria
    print("\nExtracting Bacteria...")
    bacteria_subtree = extract_domain_subtree(newick_text, "Bacteria")
    if bacteria_subtree:
        bacteria_dir = script_dir / "bacteria"
        bacteria_dir.mkdir(exist_ok=True)
        output_file = bacteria_dir / f"bacteria_tree_{timestamp}.nwk"
        with open(output_file, 'w') as f:
            f.write(bacteria_subtree)
        print(f"  Created: {output_file.name}")
        print(f"  Genera: {count_leaves(bacteria_subtree)}")
        print(f"  Size: {len(bacteria_subtree):,} characters")
    else:
        print("  ERROR: Could not extract Bacteria subtree")

if __name__ == "__main__":
    main()

