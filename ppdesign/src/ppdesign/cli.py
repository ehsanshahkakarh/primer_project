#!/usr/bin/env python3
"""
Top-level PPDesign CLI aggregator.

Provides a single entry point with subcommands:
- ppdesign unified
- ppdesign nucleotide
- ppdesign primer
- ppdesign grna
- ppdesign select
- ppdesign rank
"""

import typer

from .probedesign_unified import app as unified_app
from .probedesign_nucleotide import app as nucleotide_app
from .probedesign_primer import app as primer_app
from .probedesign_grna import app as grna_app
from .probedesign_seqselection import app as select_app
from .probedesign_rank import app as rank_app


app = typer.Typer(help="PPDesign: probe/primer and gRNA design toolkit")

# Attach existing apps as subcommands
app.add_typer(unified_app, name="unified", help="Gene-based probe/primer pipeline")
app.add_typer(nucleotide_app, name="nucleotide", help="Direct nucleotide pipeline")
app.add_typer(primer_app, name="primer", help="PCR primer pair design")
app.add_typer(grna_app, name="grna", help="CRISPR gRNA design")
app.add_typer(select_app, name="select", help="Filter and score candidate oligos")
app.add_typer(rank_app, name="rank", help="Rank probes against a reference genome")


if __name__ == "__main__":
    app()

