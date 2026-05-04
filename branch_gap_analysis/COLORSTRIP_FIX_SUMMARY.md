# Tree Colorstrip Fix Summary

**Date:** 2026-05-04  
**Issue:** Huge disconnect between tree color strips and actual novelty factor data

---

## Problem Diagnosed

### The Issue
The 18S genus tree colorstrips (division, family, genus) were showing **aggregated parent taxon data** instead of **genus-specific novelty factors**. This caused a major disconnect:

1. **Division colorstrip** was looking up the parent division (e.g., "Discosea", "Tubulinea") and showing its aggregated NF
2. **Family colorstrip** was looking up the parent family and showing its aggregated NF  
3. Only the **genus colorstrip** was showing actual genus-level data

### Specific Examples

**Before Fix:**

| Genus | Division Strip (OLD) | Genus Strip | Actual Genus NF |
|-------|---------------------|-------------|-----------------|
| Arcellinida.U.genus | 1323.0 (Tubulinea division) | inf | inf |
| Amoebozoa.U.genus | #AAAAAA (no data) | 66.9 | 66.9 |
| Acanthamoeba | 14.8 (Discosea division) | 5.0 | 5.0 |
| Amphitrema | 7.6 (Stramenopiles) | inf | inf |

### Why This Happened

The `create_itol_annotations.py` script was:
1. Loading **three separate final_merger files**:
   - `18s_ncbi_merged_division.csv` - Only 22 high-level taxa (Opisthokonta, Stramenopiles, etc.)
   - `18s_ncbi_merged_family.csv` - Only 314 family-level taxa
   - `18s_ncbi_merged_genus.csv` - All 491 genera

2. For each genus in the tree, it would:
   - **Division strip**: Walk up the lineage, find a division-rank ancestor, look it up → **WRONG** (aggregated data)
   - **Family strip**: Walk up the lineage, find a family-rank ancestor, look it up → **WRONG** (aggregated data)
   - **Genus strip**: Look up the genus directly → **CORRECT**

3. The problem: **Aggregated parent data doesn't reflect genus-specific novelty**
   - Example: Tubulinea division has 1,323 OTUs and 2 NCBI genomes (NF=1323)
   - But Arcellinida.U.genus (a genus in Tubulinea) has 75 OTUs and 0 NCBI genomes (NF=inf)
   - The division strip showed 1323, but the actual genus has NO genomic data at all!

---

## Solution Implemented

### The Fix
**All three colorstrips now use GENUS-LEVEL novelty factors** from `18s_ncbi_merged_genus.csv`.

### Code Changes

**File:** `primer_project/branch_gap_analysis/src/create_itol_annotations.py`

**Key changes:**

1. **Simplified lookup function** (lines 312-324):
```python
def _nf_for_row(level: str, genus_name: str, data: dict, nf_lookup: dict) -> float | None:
    """
    Get novelty factor for a genus at a given taxonomic level.
    
    IMPORTANT: Always returns GENUS-LEVEL NF, not aggregated parent data.
    """
    # Always look up the genus itself, not its ancestors
    return nf_lookup.get(genus_name)
```

2. **Unified data source** (lines 580-596):
```python
# Load genus-level NF data once - used for ALL colorstrips
genus_merger_file = FINAL_MERGER_FILES['genus']
genus_nf_lookup = load_nf_lookup(genus_merger_file, 'genus')

# Create all three colorstrips using genus-level data
for level in ['division', 'family', 'genus']:
    output_file = OUTPUT_DIR / f"18s_{level}_colorstrip.txt"
    create_colorstrip(level, genus_data, genus_nf_lookup, output_file)
```

3. **Updated labels** to clarify the data source:
   - Division strip: "18S Division (genus-level NF)"
   - Family strip: "18S Family (genus-level NF)"
   - Genus strip: "18S Genus"

---

## Verification

**After Fix:**

| Genus | Division Strip (NEW) | Genus Strip | Match? |
|-------|---------------------|-------------|--------|
| Arcellinida.U.genus | inf | inf | ✓ FIXED |
| Amoebozoa.U.genus | 66.9 | 66.9 | ✓ FIXED |
| Acanthamoeba | 5.0 | 5.0 | ✓ FIXED |
| Amphitrema | inf | inf | ✓ FIXED |

**Statistics:**
- All 396 tree nodes now show consistent genus-level NF across all three colorstrips
- 0 nodes show grey (no data) incorrectly
- 265 genera have NCBI genomic data (finite NF)
- 226 genera are census-only (inf NF)

---

## Impact

### What Changed
✅ **Division colorstrip** now shows **actual genus-level novelty**, not parent division aggregates  
✅ **Family colorstrip** now shows **actual genus-level novelty**, not parent family aggregates  
✅ **Genus colorstrip** unchanged (already correct)

### What This Means
- You can now trust **all three colorstrips** to show the actual novelty factor for each genus
- Colors accurately reflect whether each **specific genus** has genomic representation
- No more misleading colors where a genus appears well-represented (green/orange) but actually has no NCBI data

### Files Updated
- ✅ `primer_project/branch_gap_analysis/src/create_itol_annotations.py` - Fixed colorstrip generation
- ✅ `primer_project/branch_gap_analysis/output/18s/tree/18s_division_colorstrip.txt` - Regenerated
- ✅ `primer_project/branch_gap_analysis/output/18s/tree/18s_family_colorstrip.txt` - Regenerated
- ✅ `primer_project/branch_gap_analysis/output/18s/tree/18s_genus_colorstrip.txt` - Regenerated

---

## Tree Pruning (Optional - Not Implemented)

A pruning script was created (`prune_tree.py`) but **not executed** per your request. It can optionally:
- Remove census-only genera (inf NF)
- Remove genera with very high NF (>100)
- Create a focused tree showing only genera with genomic representation

To use it later if desired:
```bash
cd primer_project/branch_gap_analysis/src
python prune_tree.py              # Standard: remove inf NF only
python prune_tree.py --aggressive # Also remove NF > 100
```

---

## Next Steps

✅ The colorstrip data now accurately reflects genus-level novelty factors  
✅ Upload `eukcensus_18s_genus_tree.nwk` to iTOL with the updated colorstrip files  
✅ All three strips will now show consistent, accurate genus-level data  

The tree still contains all genera - nothing has been removed.
