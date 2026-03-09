# Quick Reference: Multi-Size Primer Design

## 🎯 Three Product Sizes

| Size | Range | Best For |
|------|-------|----------|
| **500bp** | 450-550bp | Degraded DNA, qPCR, NGS |
| **1000bp** | 950-1050bp | Sanger sequencing, standard PCR |
| **1500bp** | 1450-1550bp | Phylogenetics, maximum coverage |

## 🔧 Primer Settings (All Sizes)

```
Tm:     55-65°C (optimal: 60°C)
GC:     40-60%
Length: 18-25bp (optimal: 20bp)
Return: 5 primer pairs per size
```

## 📁 Output Files

```
{taxon}/primers/
├── {taxon}_consensus.fasta           # Consensus sequence
├── {taxon}_primers_all_sizes.json    # All results combined
├── {taxon}_primers_500bp.json        # 500bp primers only
├── {taxon}_primers_1000bp.json       # 1000bp primers only
└── {taxon}_primers_1500bp.json       # 1500bp primers only
```

## 💻 Usage

### Single Taxon
```bash
python scripts/design_primers.py \
  --alignment path/to/alignment.fasta \
  --output output_dir/
```

### All Taxa
```bash
python scripts/design_primers.py --gene 18S --all
```

### Test (No primer3 required)
```bash
python scripts/test_multi_size_primers.py
```

## 📊 Example Output

```json
{
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
  "1000bp": { ... },
  "1500bp": { ... }
}
```

## 🔍 Reading Results

### Python
```python
import json

# Load all sizes
with open('taxon_primers_all_sizes.json') as f:
    data = json.load(f)
    
# Get 500bp primers
primers_500 = data['500bp']['primers']
best_500 = primers_500[0]  # Best primer pair

print(f"Forward: {best_500['left_sequence']}")
print(f"Reverse: {best_500['right_sequence']}")
print(f"Product: {best_500['product_size']} bp")

# Or load specific size
with open('taxon_primers_500bp.json') as f:
    data_500 = json.load(f)
    primers = data_500['primers']
```

## ⚠️ Troubleshooting

### Sequence Too Short
If consensus < target size, that size will be skipped:
```json
{
  "1500bp": {
    "error": "Sequence too short (850 bp) for 1500bp product",
    "primers": []
  }
}
```

### No Primers Found
Check:
- Sequence has conserved regions
- Sequence is long enough
- GC content is within range (40-60%)

## 📚 Documentation

- **Full details**: `MULTI_SIZE_PRIMERS.md`
- **Changes**: `CHANGES_SUMMARY.md`
- **General usage**: `README.md`

