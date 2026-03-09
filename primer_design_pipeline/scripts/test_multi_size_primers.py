#!/usr/bin/env python3
"""
Test script to demonstrate the multi-size primer design functionality.
This shows what the new design_primers.py will do for 500bp, 1000bp, and 1500bp products.
"""

import json
from pathlib import Path

def show_multi_size_design():
    """Demonstrate the multi-size primer design approach."""
    
    print("=" * 70)
    print("Multi-Size Primer Design Configuration")
    print("=" * 70)
    print()
    
    # Define the three product sizes
    product_sizes = [500, 1000, 1500]
    
    print("Target Product Sizes:")
    print("-" * 70)
    
    for target_size in product_sizes:
        # Calculate the range (±50bp around target)
        min_size = max(target_size - 50, 100)
        max_size = target_size + 50
        
        print(f"\n{target_size}bp Product:")
        print(f"  • Target size: {target_size} bp")
        print(f"  • Allowed range: {min_size}-{max_size} bp")
        print(f"  • Primer3 will search for primers that produce amplicons in this range")
        print(f"  • Returns: Top 5 primer pairs for this size")
    
    print("\n" + "=" * 70)
    print("Primer3 Settings (Same for All Sizes)")
    print("=" * 70)
    
    settings = {
        "Primer Length": "18-25 bp (optimal: 20 bp)",
        "Melting Temperature (Tm)": "55-65°C (optimal: 60°C)",
        "GC Content": "40-60%",
        "Number of Pairs Returned": "5 per product size"
    }
    
    for param, value in settings.items():
        print(f"  • {param}: {value}")
    
    print("\n" + "=" * 70)
    print("Output Files Generated")
    print("=" * 70)
    print()
    print("For each taxon, the following files will be created:")
    print("  1. {taxon}_primers_all_sizes.json - All results combined")
    print("  2. {taxon}_primers_500bp.json - Only 500bp primers")
    print("  3. {taxon}_primers_1000bp.json - Only 1000bp primers")
    print("  4. {taxon}_primers_1500bp.json - Only 1500bp primers")
    print("  5. {taxon}_consensus.fasta - Consensus sequence used")
    
    print("\n" + "=" * 70)
    print("Example Output Structure")
    print("=" * 70)
    
    example_output = {
        "500bp": {
            "target_size": 500,
            "size_range": [450, 550],
            "primers": [
                {
                    "pair_id": 0,
                    "left_sequence": "CTCATCCGAAACGCTGCATG",
                    "right_sequence": "TTGGCTGCGGATTACCCATT",
                    "left_tm": 60.2,
                    "right_tm": 59.8,
                    "left_gc": 50.0,
                    "right_gc": 55.0,
                    "product_size": 495,
                    "penalty": 0.234
                }
            ]
        },
        "1000bp": {
            "target_size": 1000,
            "size_range": [950, 1050],
            "primers": [
                {
                    "pair_id": 0,
                    "product_size": 1005
                }
            ]
        },
        "1500bp": {
            "target_size": 1500,
            "size_range": [1450, 1550],
            "primers": [
                {
                    "pair_id": 0,
                    "product_size": 1498
                }
            ]
        }
    }
    
    print("\nJSON structure (simplified):")
    print(json.dumps(example_output, indent=2))
    
    print("\n" + "=" * 70)
    print("Usage")
    print("=" * 70)
    print()
    print("To design primers for all three sizes:")
    print("  python scripts/design_primers.py --alignment path/to/alignment.fasta")
    print()
    print("To process all alignments in test_data:")
    print("  python scripts/design_primers.py --gene 18S --all")
    print()

if __name__ == '__main__':
    show_multi_size_design()

