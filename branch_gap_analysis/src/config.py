"""
Configuration for branch_gap_analysis pipeline.
Points to source NCBI CSV files from the OTU assembly pipeline.
"""

from pathlib import Path

# Base paths
PROJECT_ROOT = Path(__file__).parent.parent
REPO_ROOT = PROJECT_ROOT.parent
PROJA_ROOT = REPO_ROOT.parent  # /Users/ehsankakarh/PROJA

# Source data paths (NCBI merged files from GitHub repo)
OTU_PIPELINE_ROOT = PROJA_ROOT / "GitHub" / "otu_assembly_comparative_pipeline-"
SOURCE_DATA_BASE = OTU_PIPELINE_ROOT / "final_merger" / "outputs"

# 16S NCBI files
SOURCE_16S_FILES = {
    "division": SOURCE_DATA_BASE / "16s_ncbi_merged_division.csv",
    "family": SOURCE_DATA_BASE / "16s_ncbi_merged_family.csv",
    "genus": SOURCE_DATA_BASE / "16s_ncbi_merged_genus.csv",
}

# 18S NCBI files
SOURCE_18S_FILES = {
    "division": SOURCE_DATA_BASE / "18s_ncbi_merged_division.csv",
    "family": SOURCE_DATA_BASE / "18s_ncbi_merged_family.csv",
    "genus": SOURCE_DATA_BASE / "18s_ncbi_merged_genus.csv",
}

# Output paths (all outputs go to output/ directory)
OUTPUT_DIR = PROJECT_ROOT / "output"
METADATA_DIR = OUTPUT_DIR  # Legacy alias - now points to output/

# Taxonomic ranks in order (standard NCBI ranks)
TAXONOMIC_RANKS = [
    "superkingdom", "domain", "kingdom", "subkingdom",
    "superphylum", "phylum", "subphylum",
    "superclass", "class", "subclass", "infraclass",
    "superorder", "order", "suborder", "infraorder", "parvorder",
    "superfamily", "family", "subfamily",
    "tribe", "subtribe",
    "genus", "subgenus",
    "species_group", "species_subgroup", "species",
    "subspecies", "strain"
]

# Simplified ranks for output columns
MAIN_RANKS = [
    "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"
]

# Analysis parameters
TOP_N = 10  # Number of top novelty_factor entries to analyze

