import os
import typer
import polars as pl
from Bio import AlignIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import logging
import numpy as np
from typing import Tuple
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

app = typer.Typer()
# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def process_alignment_file(
    file, codon_dir, length_range, window_size, step, threshold, gc_range, tm_range
):
    """Process a single alignment file, adapted for length_range."""
    alignment_path = os.path.join(codon_dir, file)
    logging.info(f"Reading alignment file: {file}")
    alignment = AlignIO.read(alignment_path, "fasta")
    min_len, max_len = length_range  # Unpack the length range
    return find_conserved_regions(
        alignment,
        min_len,
        max_len,
        window_size,
        step,
        threshold,
        file,
        gc_range,
        tm_range,
    )


def generate_consensus_sequence(region):
    """Generate consensus sequence with degenerate nucleotides for non-identical positions."""
    consensus = []
    for i in range(len(region[0].seq)):
        bases = [record.seq[i] for record in region]
        base_count = Counter(bases)
        if len(base_count) == 1:
            # If all bases are identical, add the base to the consensus
            consensus.append(bases[0])
        else:
            # Otherwise, add a degenerate base according to the observed bases
            degenerate_bases = {
                frozenset(["A", "G"]): "R",  # Purine
                frozenset(["C", "T"]): "Y",  # Pyrimidine
                frozenset(["G", "C"]): "S",
                frozenset(["A", "T"]): "W",
                frozenset(["G", "T"]): "K",
                frozenset(["A", "C"]): "M",
                frozenset(["C", "G", "T"]): "B",
                frozenset(["A", "G", "T"]): "D",
                frozenset(["A", "C", "T"]): "H",
                frozenset(["A", "C", "G"]): "V",
                frozenset(["A", "C", "G", "T"]): "N",  # Any base
            }
            consensus_base = degenerate_bases.get(frozenset(base_count.keys()), "N")
            consensus.append(consensus_base)
    return "".join(consensus)


def calculate_nd_perc(sequence):
    """Calculate the percentage of non-degenerate (A, T, C, G) bases in a sequence."""
    non_degenerate_bases = "ATCG"
    non_degenerate_count = sum(1 for base in sequence if base in non_degenerate_bases)
    total_bases = len(sequence)
    nd_perc = (non_degenerate_count / total_bases) * 100 if total_bases > 0 else 0
    return nd_perc


def calculate_composite_score(tm, gc_content, conservation, tm_range, gc_range):
    # Normalize each parameter to a 0-1 range based on expected ranges
    tm_normalized = (tm - tm_range[0]) / (tm_range[1] - tm_range[0])
    gc_normalized = (gc_content - gc_range[0]) / (gc_range[1] - gc_range[0])
    conservation_normalized = (
        conservation / 100
    )  # Assuming conservation is already a percentage

    # Calculate composite score as an average of normalized values
    composite_score = (tm_normalized + gc_normalized + conservation_normalized) / 3
    return composite_score


def parse_range(ctx, param, value):
    """Directly returns the tuple without parsing a string."""
    # Since Typer now correctly handles tuple input, we directly return the value.
    return value


def calculate_tm_gc(sequence: str) -> Tuple[float, float]:
    """Calculates the melting temperature (Tm) and GC content for a given DNA sequence."""
    seq = Seq(sequence)
    tm = mt.Tm_Wallace(seq, strict=False)  # Wallace's rule: 4*(G+C) + 2*(A+T) in °C
    gc_content = (sequence.count("G") + sequence.count("C")) / len(sequence) * 100
    return tm, gc_content


def is_within_range(value: float, value_range: Tuple[float, float]) -> bool:
    """Checks if a value is within a specified range."""
    return value_range[0] <= value <= value_range[1]


def reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence, ignoring non-ATCG characters."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    reversed_complement = []
    for base in reversed(seq):
        if base in complement:
            reversed_complement.append(complement[base])
        else:
            reversed_complement.append("N")
            # continue
    return "".join(reversed_complement)


def check_hairpin(sequence, min_loop=3, min_stem=4):
    """Checks for hairpin structures in a sequence."""
    for i in range(len(sequence) - 2 * min_stem - min_loop):
        stem = sequence[i : i + min_stem]
        loop_start = i + min_stem + min_loop
        for j in range(loop_start, len(sequence) - min_stem + 1):
            if sequence[j : j + min_stem] == reverse_complement(stem):
                return True  # Hairpin structure detected
    return False


def check_dimer(sequence, min_dimer_len=4):
    """Checks for dimer formation within a sequence."""
    rev_comp_seq = reverse_complement(sequence)
    for i in range(len(sequence) - min_dimer_len + 1):
        for j in range(len(rev_comp_seq) - min_dimer_len + 1):
            if sequence[i : i + min_dimer_len] == rev_comp_seq[j : j + min_dimer_len]:
                return True  # Potential dimer formation detected
    return False


def calculate_positional_conservation(region):
    """Calculate conservation scores, penalizing for nucleotide diversity."""
    num_sequences = len(region)
    conservation_scores = []
    penalty_factor = 0.05

    for i in range(len(region[0].seq)):
        nucleotides = [record.seq[i] for record in region]
        nucleotide_count = Counter(nucleotides)
        most_common_nucleotide, most_common_count = nucleotide_count.most_common(1)[0]
        # Adjust conservation score based on diversity
        # Example: decrease score based on the number of different nucleotides present
        diversity_penalty = (
            len(nucleotide_count) - 1
        )  # Number of non-most-common nucleotides
        adjusted_score = (most_common_count / num_sequences) - (
            diversity_penalty * penalty_factor
        )
        conservation_scores.append(
            max(0, adjusted_score)
        )  # Ensure score doesn't go negative

    return conservation_scores


def find_conserved_regions(
    alignment,
    min_len,
    max_len,
    window_size,
    step,
    threshold,
    file_name,
    gc_range,
    tm_range,
):
    conserved_regions = []
    for window_start in range(0, len(alignment[0]) - window_size + 1, step):
        for region_start in range(
            window_start, window_start + window_size - max_len + 1
        ):
            for region_end in range(
                region_start + min_len,
                min(region_start + max_len + 1, window_start + window_size),
            ):
                region = alignment[:, region_start:region_end]
                conservation_scores = calculate_positional_conservation(region)
                if all(score >= threshold for score in conservation_scores):
                    sequence = generate_consensus_sequence(region)
                    tm, gc_content = calculate_tm_gc(sequence)
                    hairpin_check = check_hairpin(sequence)
                    dimer_check = check_dimer(sequence)
                    if (
                        not hairpin_check
                        and not dimer_check
                        and is_within_range(gc_content, gc_range)
                        and is_within_range(tm, tm_range)
                    ):
                        gene_id = (
                            region[0].id.split("|")[1]
                            if "|" in region[0].id
                            else region[0].id
                        )
                        conservation = np.mean(conservation_scores) * 100
                        composite_score = calculate_composite_score(
                            tm, gc_content, conservation, tm_range, gc_range
                        )
                        nd_perc = calculate_nd_perc(sequence)
                        conserved_regions.append(
                            {
                                "og_id": file_name.replace(".codon", ""),
                                "gene_id": gene_id,
                                "start": region_start,
                                "end": region_end,
                                "length": region_end - region_start,
                                "sequence": sequence,
                                "tm": tm,
                                "gc%": gc_content,
                                "conservation": conservation,
                                "composite_score": composite_score,
                                "nd_perc": nd_perc,  # Add the nd_perc value here
                            }
                        )
    return conserved_regions


@app.command()
def main(
    codon_dir: str = typer.Option(
        ..., "-i", "--input", help="Directory containing codon alignments"
    ),
    length_range: Tuple[int, int] = typer.Option(
        (10, 12),
        "-l",
        "--length",
        help="Range of lengths for conserved regions",
        callback=parse_range,
    ),
    window_size: int = typer.Option(20, "--window-size", help="Sliding window size"),
    step: int = typer.Option(5, "--step", help="Step size for sliding window"),
    threshold: float = typer.Option(
        0.9, "--threshold", help="Threshold for conservation"
    ),
    gc_range: Tuple[float, float] = typer.Option(
        (40.0, 60.0), "--gc-range", help="GC content range (%)", callback=parse_range
    ),
    tm_range: Tuple[float, float] = typer.Option(
        (50.0, 60.0),
        "--tm-range",
        help="Melting temperature (Tm) range (°C)",
        callback=parse_range,
    ),
    threads: int = typer.Option(
        4, "--threads", help="Number of threads for parallel processing"
    ),
):
    conserved_regions_info = (
        []
    )  # This will store dictionaries for each conserved region found

    logging.info(
        f"Processing codon alignments in directory: {codon_dir} using {threads} threads"
    )
    logging.info(f"Using GC content range: {gc_range}% and Tm range: {tm_range}°C")
    logging.info(f"Filtering for sequence lengths between: {length_range}")

    # Use ThreadPoolExecutor to process files in parallel
    files = [file for file in os.listdir(codon_dir) if file.endswith(".codon")]
    # Initialize tqdm with the total number of files to process
    with (
        ThreadPoolExecutor(max_workers=threads) as executor,
        tqdm(total=len(files), desc="Processing Files") as progress,
    ):
        future_to_file = {
            executor.submit(
                process_alignment_file,
                file,
                codon_dir,
                length_range,
                window_size,
                step,
                threshold,
                gc_range,
                tm_range,
            ): file
            for file in files
        }
        for future in as_completed(future_to_file):
            try:
                conserved_regions = future.result()
                conserved_regions_info.extend(conserved_regions)
            except Exception as exc:
                file = future_to_file[future]
                logging.error(f"{file} generated an exception: {exc}")
            finally:
                # Update the progress bar upon completion of each task
                progress.update(1)

    if conserved_regions_info:
        lengths = [len(region["sequence"]) for region in conserved_regions_info]
        logging.info(f"Minimum sequence length found: {min(lengths)}")
        logging.info(f"Maximum sequence length found: {max(lengths)}")
    else:
        logging.warning("No conserved regions found. Exiting script.")
        return

    conserved_regions_info.sort(key=lambda x: x["nd_perc"], reverse=True)
    # Convert to DataFrame
    df = pl.DataFrame(conserved_regions_info).unique()
    if df.is_empty():
        logging.warning("DataFrame is empty after creation. No data to process.")
        return  # Similar early return if the DataFrame happens to be empty
    logging.info(f"DataFrame columns: {df.columns}")
    # Before filtering, optionally log the DataFrame's structure
    logging.info(f"DataFrame structure before filtering: {df.schema}")
    df_filtered = df.filter(
        (df["gc%"] >= gc_range[0])
        & (df["gc%"] <= gc_range[1])
        & (df["tm"] >= tm_range[0])
        & (df["tm"] <= tm_range[1])
    )

    # Ensure the DataFrame is being correctly populated with data
    if df.shape[0] > 0:
        logging.info("DataFrame is populated with conserved region data.")
    else:
        logging.warning("DataFrame is empty. No conserved regions were found or added.")

    try:
        if not df.is_empty():
            df_filtered = df.filter(
                (df["gc%"] >= gc_range[0])
                & (df["gc%"] <= gc_range[1])
                & (df["tm"] >= tm_range[0])
                & (df["tm"] <= tm_range[1])
            )
            if not df_filtered.is_empty():
                df_filtered.write_csv("filtered_conserved_regions.csv")
                logging.info("Filtered conserved regions extracted and saved.")
            else:
                logging.warning(
                    "No conserved regions found within the specified GC% and Tm ranges."
                )
        else:
            logging.warning("DataFrame is empty. Skipping filtering step.")
    except Exception as e:
        logging.error(f"Error filtering DataFrame: {e}")


if __name__ == "__main__":
    app()
