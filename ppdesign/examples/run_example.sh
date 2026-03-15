#!/bin/bash
# Example workflow for PPDesign

# Set up environment
echo "Setting up PPDesign environment..."
pixi shell

# Example 1: Design probes for a viral family using GeNomad output
echo "Example 1: Designing probes for Faserviricetes..."
ppdesign-genomad \
  -t ../data/test/genomad_virus_summary.tsv \
  -f ../data/test/genomad_virus.fna \
  -x "Faserviricetes" \
  -c 0.3 \
  -j 4 \
  -r \
  -n 10  # Using small number for example

# Example 2: Select oligonucleotides
echo "Example 2: Selecting oligonucleotides..."
ppdesign-select \
  -i Faserviricetes_og30_codon \
  --gc-range 30 60 \
  --tm-range 35 45 \
  --length 15 25 \
  --window-size 25 \
  --step 5 \
  --threshold 0.8 \
  --threads 4

# Example 3: Rank and visualize probes
echo "Example 3: Ranking probes..."
ppdesign-rank \
  -o Faserviricetes_og30_codon_selected_oligos.csv \
  -g ../data/test/reference_genome.fna \
  -r probe_ranking_results

echo "Example workflow completed!"
echo "Check the output files for results."