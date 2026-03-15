"""PPDesign: Probe and Primer Design for Targeted Taxonomic Groups.

A bioinformatics pipeline for designing oligonucleotide probes and primers
targeting specific taxonomic groups with support for degenerate nucleotides
and comprehensive thermodynamic filtering.
"""

__version__ = "0.9.0"
__author__ = "PPDesign Development Team"

from . import (
    conserved_finder,
    genome_database,
    kmer_finder,
    probedesign_nucleotide,
    probedesign_rank,
    probedesign_seqselection,
    probedesign_unified,
)

__all__ = [
    "conserved_finder",
    "genome_database",
    "kmer_finder",
    "probedesign_nucleotide",
    "probedesign_rank",
    "probedesign_seqselection",
    "probedesign_unified",
]
