import pandas as pd
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
import argparse


def canonical_kmer(kmer):
    """Return the lexicographically smallest k-mer among the k-mer and its reverse complement."""
    seq = Seq(kmer)
    rev_comp_kmer = str(seq.reverse_complement())
    # Only consider forward and reverse complement (standard practice for canonical k-mers)
    variants = [kmer, rev_comp_kmer]
    canonical = min(variants)
    return canonical


def count_kmers(sequence, k_values):
    counts_per_k = {k: Counter() for k in k_values}
    seq_len = len(sequence)
    for k in k_values:
        counts = counts_per_k[k]
        for i in range(seq_len - k + 1):
            kmer = str(sequence[i : i + k]).upper()
            if "N" in kmer:
                continue  # Skip kmers containing ambiguous bases
            canonical = canonical_kmer(kmer)
            counts[canonical] += 1
    return counts_per_k


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Count kmers in a FASTA file considering reverse and reverse complement"
    )
    parser.add_argument(
        "fasta_file", help="Input FASTA file containing nucleotide sequences"
    )
    parser.add_argument("output_csv", help="Output CSV file to store k-mer counts")
    parser.add_argument(
        "--kmer-lengths",
        default="4,5",
        help="Comma-separated k-mer lengths (e.g., 4,5,6)",
    )
    args = parser.parse_args()

    # Parse k-mer lengths
    k_values = [int(k) for k in args.kmer_lengths.split(",")]

    contig_kmer_counts = {}
    total_kmer_counts_per_k = {k: Counter() for k in k_values}
    contig_sequence_lengths = {}

    # Parse the FASTA file using Biopython
    for record in SeqIO.parse(args.fasta_file, "fasta"):
        contig_id = record.id
        sequence = record.seq
        seq_len = len(sequence)
        counts_per_k = count_kmers(sequence, k_values)
        contig_kmer_counts[contig_id] = counts_per_k
        contig_sequence_lengths[contig_id] = seq_len
        for k in k_values:
            total_kmer_counts_per_k[k].update(counts_per_k[k])

    # Identify the 10 most frequent kmers for each k
    top_kmers_per_k = {}
    for k in k_values:
        top_kmers = [kmer for kmer, _ in total_kmer_counts_per_k[k].most_common(10)]
        top_kmers_per_k[k] = top_kmers

    # Prepare columns for DataFrame
    columns = []
    for k in k_values:
        for kmer in top_kmers_per_k[k]:
            columns.append(kmer)
            columns.append(f"{kmer}_perc")

    # Prepare data for DataFrame
    data = []
    contig_ids = []
    for contig_id in contig_sequence_lengths:
        seq_len = contig_sequence_lengths[contig_id]
        counts_per_k = contig_kmer_counts[contig_id]
        row = []
        for k in k_values:
            counts = counts_per_k[k]
            for kmer in top_kmers_per_k[k]:
                count = counts.get(kmer, 0)
                percentage = (count * k) / seq_len * 100 if seq_len > 0 else 0
                row.append(count)
                row.append(percentage)
        data.append(row)
        contig_ids.append(contig_id)

    # Add total counts
    total_row = []
    total_seq_len = sum(contig_sequence_lengths.values())
    for k in k_values:
        total_counts = total_kmer_counts_per_k[k]
        for kmer in top_kmers_per_k[k]:
            count = total_counts.get(kmer, 0)
            percentage = (count * k) / total_seq_len * 100 if total_seq_len > 0 else 0
            total_row.append(count)
            total_row.append(percentage)
    data.append(total_row)
    contig_ids.append("total")

    df = pd.DataFrame(data, index=contig_ids, columns=columns)
    df.to_csv(args.output_csv)


if __name__ == "__main__":
    main()
