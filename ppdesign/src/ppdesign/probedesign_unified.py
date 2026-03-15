#!/usr/bin/env python3
"""
Unified probe design pipeline for bacteria and viruses using gene calling.
Replaces both probedesign_codonaln.py and probedesign_codonaln_genomad.py
"""

import shutil
import logging
import subprocess
from pathlib import Path
from typing import List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import typer
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .genome_database import GenomeDatabase
from .output_utils import resolve_output_dir

app = typer.Typer()

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler("pipeline.log"), logging.StreamHandler()],
)


class UnifiedProbeDesign:
    """Unified pipeline for probe design with gene calling."""

    def _check_dependencies(self):
        """Ensure required external tools are available."""
        required = ["prodigal-gv", "proteinortho", "mafft"]
        missing = [tool for tool in required if shutil.which(tool) is None]
        if missing:
            raise RuntimeError(
                "Missing required tools for unified pipeline: {}".format(', '.join(missing))
            )

    def __init__(self, output_dir: str, threads: int = 4):
        """Initialize the pipeline."""
        self.output_dir = resolve_output_dir(output_dir)
        self.threads = threads
        self.logger = logging.getLogger(__name__)

        self._check_dependencies()

        # Create output directory structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.fna_dir = self.output_dir / "fna"
        self.faa_dir = self.output_dir / "faa"
        self.fnn_dir = self.output_dir / "fnn"
        self.proteinortho_dir = self.output_dir / "proteinortho"
        self.conserved_faa_dir = self.output_dir / "conserved_faa"
        self.conserved_fnn_dir = self.output_dir / "conserved_fnn"
        self.alignment_dir = self.output_dir / "alignments"
        self.codon_dir = self.output_dir / "codon_alignments"

    def run_database_mode(
        self,
        taxonomy: str,
        tax_level: str = "species",
        limit: Optional[int] = None,
        conservation_level: float = 0.3,
    ):
        """Run pipeline using database query mode."""
        self.logger.info(f"Running database mode for taxonomy: {taxonomy}")

        # Create data subdirectory for downloaded genomes
        data_dir = self.output_dir / "data"
        data_dir.mkdir(exist_ok=True)

        # Query and download genomes
        db = GenomeDatabase()

        # Support both bacteria and viruses
        genomes_df = db.query_by_taxonomy(
            taxonomy=taxonomy, tax_level=tax_level, limit=limit
        )

        if genomes_df.empty:
            raise ValueError(f"No genomes found for taxonomy: {taxonomy}")

        self.logger.info(f"Found {len(genomes_df)} genomes to download")

        # Download genomes (only FNA files, we'll generate FAA/FNN with prodigal-gv)
        downloaded_files = db.download_genomes(
            genomes_df=genomes_df,
            output_dir=str(data_dir),
            download_fna=True,
            download_faa=False,  # Don't download FAA, we'll generate with prodigal-gv
            threads=self.threads,
        )

        # Move files to standard directories
        self.fna_dir.mkdir(exist_ok=True)
        self.faa_dir.mkdir(exist_ok=True)

        for fna_file in downloaded_files["fna"]:
            shutil.copy2(fna_file, self.fna_dir / Path(fna_file).name)

        # Since we didn't download FAA files, we need to generate them with prodigal-gv
        self._run_gene_calling()

        # Continue with standard pipeline
        self._run_pipeline(conservation_level)

    def run_local_mode(self, fasta_input: str, conservation_level: float = 0.3):
        """Run pipeline using local FASTA file(s). Accepts a file or directory."""
        self.logger.info(f"Running local mode with input: {fasta_input}")

        fasta_path = Path(fasta_input)
        if not fasta_path.exists():
            raise ValueError(f"FASTA path not found: {fasta_input}")

        # Copy FASTA files to fna directory
        self.fna_dir.mkdir(exist_ok=True)

        if fasta_path.is_file():
            fasta_files = [fasta_path]
        else:
            fasta_files = (
                list(fasta_path.glob("*.fna"))
                + list(fasta_path.glob("*.fasta"))
                + list(fasta_path.glob("*.fa"))
            )

        if not fasta_files:
            raise ValueError(f"No FASTA files found in: {fasta_input}")

        self.logger.info(f"Found {len(fasta_files)} FASTA files")

        for fasta_file in fasta_files:
            # Standardize to .fna extension
            dest_name = fasta_file.stem + ".fna"
            shutil.copy2(fasta_file, self.fna_dir / dest_name)

        # Run gene calling
        self._run_gene_calling()

        # Continue with standard pipeline
        self._run_pipeline(conservation_level)

    def _run_gene_calling(self):
        """Run prodigal-gv on all FNA files to generate FAA and FNN."""
        self.logger.info("Running gene calling with prodigal-gv")

        self.faa_dir.mkdir(exist_ok=True)
        self.fnn_dir.mkdir(exist_ok=True)

        fna_files = list(self.fna_dir.glob("*.fna"))

        def run_prodigal_single(fna_file: Path) -> bool:
            """Run prodigal-gv on a single file."""
            base_name = fna_file.stem
            faa_output = self.faa_dir / f"{base_name}.faa"
            fnn_output = self.fnn_dir / f"{base_name}.fnn"

            cmd = [
                "prodigal-gv",
                "-i",
                str(fna_file),
                "-a",
                str(faa_output),
                "-d",
                str(fnn_output),
                "-p",
                "meta",
                "-q",
            ]

            try:
                subprocess.run(cmd, check=True, capture_output=True)

                # Adjust headers to include genome ID
                self._adjust_prodigal_headers(faa_output, fnn_output, base_name)

                return True
            except subprocess.CalledProcessError as e:
                self.logger.error(
                    f"Prodigal-gv failed on {fna_file}: {e.stderr.decode()}"
                )
                return False

        # Run in parallel
        success_count = 0
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            futures = {
                executor.submit(run_prodigal_single, fna_file): fna_file
                for fna_file in fna_files
            }

            for future in tqdm(
                as_completed(futures), total=len(futures), desc="Running prodigal-gv"
            ):
                if future.result():
                    success_count += 1

        self.logger.info(
            f"Gene calling completed: {success_count}/{len(fna_files)} successful"
        )

    def _generate_fnn_files(self):
        """Generate FNN files from existing FAA files using prodigal-gv."""
        self.logger.info("Generating FNN files")

        self.fnn_dir.mkdir(exist_ok=True)

        fna_files = list(self.fna_dir.glob("*.fna"))

        def generate_fnn_single(fna_file: Path) -> bool:
            base_name = fna_file.stem
            fnn_output = self.fnn_dir / f"{base_name}.fnn"

            if fnn_output.exists():
                return True

            cmd = [
                "prodigal-gv",
                "-i",
                str(fna_file),
                "-d",
                str(fnn_output),
                "-p",
                "meta",
                "-q",
            ]

            try:
                subprocess.run(cmd, check=True, capture_output=True)

                # Adjust headers
                self._adjust_fnn_headers(fnn_output, base_name)

                return True
            except subprocess.CalledProcessError as e:
                self.logger.error(
                    f"Failed to generate FNN for {fna_file}: {e.stderr.decode()}"
                )
                return False

        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            futures = {
                executor.submit(generate_fnn_single, fna_file): fna_file
                for fna_file in fna_files
            }

            for future in tqdm(
                as_completed(futures), total=len(futures), desc="Generating FNN files"
            ):
                future.result()

    def _adjust_prodigal_headers(self, faa_file: Path, fnn_file: Path, genome_id: str):
        """Adjust prodigal output headers to include genome ID and remove stop codons."""
        # Adjust FAA
        records = []
        for record in SeqIO.parse(faa_file, "fasta"):
            record.id = f"{genome_id}|{record.id}"
            record.description = ""
            # Remove stop codons (*) from protein sequences
            record.seq = record.seq.rstrip("*")
            records.append(record)
        SeqIO.write(records, faa_file, "fasta")

        # Adjust FNN
        records = []
        for record in SeqIO.parse(fnn_file, "fasta"):
            record.id = f"{genome_id}|{record.id}"
            record.description = ""
            records.append(record)
        SeqIO.write(records, fnn_file, "fasta")

    def _adjust_fnn_headers(self, fnn_file: Path, genome_id: str):
        """Adjust FNN headers to include genome ID."""
        records = []
        for record in SeqIO.parse(fnn_file, "fasta"):
            if "|" not in record.id:
                record.id = f"{genome_id}|{record.id}"
                record.description = ""
            records.append(record)
        SeqIO.write(records, fnn_file, "fasta")

    def _run_pipeline(self, conservation_level: float):
        """Run the main pipeline after gene calling."""
        # Step 1: Run ProteinOrtho
        self._run_proteinortho()

        # Step 2: Parse ProteinOrtho output
        proteinortho_tsv = (
            self.proteinortho_dir / f"{self.output_dir.name}.proteinortho.tsv"
        )
        orthogroups_file = (
            self.proteinortho_dir / f"{self.output_dir.name}.proteinortho.groups"
        )

        if not proteinortho_tsv.exists():
            raise FileNotFoundError(
                f"ProteinOrtho output not found: {proteinortho_tsv}"
            )

        self._create_orthogroups_file(proteinortho_tsv, orthogroups_file)

        # Step 3: Filter orthogroups by conservation
        total_genomes = len(list(self.faa_dir.glob("*.faa")))
        qualifying_orthogroups = self._filter_orthogroups_by_conservation(
            orthogroups_file, conservation_level, total_genomes
        )

        self.logger.info(
            f"Found {len(qualifying_orthogroups)} orthogroups meeting conservation level {conservation_level}"
        )

        if not qualifying_orthogroups:
            self.logger.warning("No conserved orthogroups found")
            return

        # Step 4: Create multi-FAA and FNN files
        self._create_multi_faa_fnn_files(qualifying_orthogroups, orthogroups_file)

        # Step 5: Align with MAFFT
        self._align_all_faa_files()

        # Step 6: Create codon alignments
        self._create_codon_alignments(qualifying_orthogroups)

        self.logger.info("Pipeline completed successfully")

    def _run_proteinortho(self):
        """Run ProteinOrtho on all FAA files."""
        self.logger.info("Running ProteinOrtho")

        self.proteinortho_dir.mkdir(exist_ok=True)

        faa_files = list(self.faa_dir.glob("*.faa"))

        if not faa_files:
            raise ValueError("No FAA files found for ProteinOrtho")

        # Run proteinortho with absolute paths
        project_name = self.output_dir.name
        cmd = ["proteinortho", f"-project={project_name}", f"-cpus={self.threads}"] + [
            str(f.absolute()) for f in faa_files
        ]

        try:
            # Run in the proteinortho directory
            subprocess.run(
                cmd, check=True, capture_output=True, cwd=str(self.proteinortho_dir)
            )
            self.logger.info("ProteinOrtho completed successfully")

        except subprocess.CalledProcessError as e:
            self.logger.error(f"ProteinOrtho failed: {e.stderr.decode()}")
            raise

    def _create_orthogroups_file(self, proteinortho_tsv: Path, output_file: Path):
        """Parse ProteinOrtho output to create orthogroups file."""
        with open(proteinortho_tsv, "r") as infile, open(output_file, "w") as outfile:
            header = infile.readline().strip().split("\t")[3:]  # Skip initial columns

            for i, line in enumerate(infile):
                parts = line.strip().split("\t")
                orthogroup = f"og{i+1}"
                proteins = []

                for genome, protein_list in zip(header, parts[3:]):
                    genome_id = Path(genome).stem

                    if protein_list and not protein_list.startswith("*"):
                        for protein in protein_list.split(","):
                            protein_id = (
                                protein.split("|")[-1] if "|" in protein else protein
                            )
                            proteins.append(f"{genome_id}|{protein_id}")

                if proteins:
                    outfile.write(f"{orthogroup} {' '.join(proteins)}\n")

    def _filter_orthogroups_by_conservation(
        self, orthogroups_file: Path, conservation_level: float, total_genomes: int
    ) -> List[str]:
        """Filter orthogroups by conservation level."""
        qualifying_orthogroups = []

        with open(orthogroups_file, "r") as infile:
            for line in infile:
                parts = line.strip().split()
                orthogroup = parts[0]

                # Count unique genomes
                unique_genomes = set(part.split("|")[0] for part in parts[1:])

                if len(unique_genomes) >= round(conservation_level * total_genomes):
                    qualifying_orthogroups.append(orthogroup)

        return qualifying_orthogroups

    def _create_multi_faa_fnn_files(
        self, qualifying_orthogroups: List[str], orthogroups_file: Path
    ):
        """Create multi-sequence FAA and FNN files for conserved orthogroups."""
        self.logger.info("Creating multi-FAA and FNN files")

        self.conserved_faa_dir.mkdir(exist_ok=True)
        self.conserved_fnn_dir.mkdir(exist_ok=True)

        # Load orthogroup data
        with open(orthogroups_file, "r") as f:
            orthogroup_data = {line.split()[0]: line.strip().split()[1:] for line in f}

        def process_orthogroup(og_id: str):
            if og_id not in orthogroup_data:
                return

            proteins = orthogroup_data[og_id]
            faa_records = []
            fnn_records = []

            for protein_info in proteins:
                genome_id, protein_id = protein_info.split("|", 1)

                # Read FAA
                faa_file = self.faa_dir / f"{genome_id}.faa"
                if faa_file.exists():
                    for record in SeqIO.parse(faa_file, "fasta"):
                        if record.id == protein_info or record.id.endswith(
                            f"|{protein_id}"
                        ):
                            faa_records.append(record)
                            break

                # Read FNN
                fnn_file = self.fnn_dir / f"{genome_id}.fnn"
                if fnn_file.exists():
                    for record in SeqIO.parse(fnn_file, "fasta"):
                        if protein_id in record.id:
                            fnn_records.append(record)
                            break

            # Write files
            if faa_records:
                SeqIO.write(
                    faa_records, self.conserved_faa_dir / f"{og_id}.faa", "fasta"
                )

            if fnn_records:
                SeqIO.write(
                    fnn_records, self.conserved_fnn_dir / f"{og_id}.fnn", "fasta"
                )

        # Process in parallel
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            futures = [
                executor.submit(process_orthogroup, og_id)
                for og_id in qualifying_orthogroups
            ]

            for future in tqdm(
                as_completed(futures),
                total=len(futures),
                desc="Creating multi-sequence files",
            ):
                future.result()

    def _align_all_faa_files(self):
        """Align all FAA files using MAFFT."""
        self.logger.info("Running MAFFT alignments")

        self.alignment_dir.mkdir(exist_ok=True)

        faa_files = list(self.conserved_faa_dir.glob("*.faa"))

        def align_single_file(faa_file: Path):
            output_file = self.alignment_dir / f"{faa_file.stem}.faa.mafft"

            cmd = ["mafft", "--auto", str(faa_file)]

            try:
                with open(output_file, "w") as out:
                    subprocess.run(
                        cmd, stdout=out, stderr=subprocess.DEVNULL, check=True
                    )
            except subprocess.CalledProcessError:
                self.logger.error(f"MAFFT failed on {faa_file}")

        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            futures = [executor.submit(align_single_file, f) for f in faa_files]

            for future in tqdm(
                as_completed(futures), total=len(futures), desc="Running alignments"
            ):
                future.result()

    def _create_codon_alignments(self, qualifying_orthogroups: List[str]):
        """Create codon alignments from protein alignments and nucleotide sequences."""
        self.logger.info("Creating codon alignments")

        self.codon_dir.mkdir(exist_ok=True)

        def create_single_codon_alignment(og_id: str):
            protein_alignment_file = self.alignment_dir / f"{og_id}.faa.mafft"
            nucleotide_file = self.conserved_fnn_dir / f"{og_id}.fnn"
            output_file = self.codon_dir / f"{og_id}.codon"

            if not protein_alignment_file.exists() or not nucleotide_file.exists():
                return

            try:
                # Load protein alignment
                protein_alignment = AlignIO.read(protein_alignment_file, "fasta")

                # Load nucleotide sequences
                nuc_seqs = {
                    record.id: str(record.seq)
                    for record in SeqIO.parse(nucleotide_file, "fasta")
                }

                # Create codon alignment
                codon_records = []

                for protein_record in protein_alignment:
                    protein_seq = str(protein_record.seq)

                    # Find matching nucleotide sequence
                    nuc_seq = None
                    for nuc_id, seq in nuc_seqs.items():
                        if protein_record.id in nuc_id or nuc_id in protein_record.id:
                            nuc_seq = seq
                            break

                    if not nuc_seq:
                        continue

                    # Build codon sequence
                    codon_seq = []
                    nuc_index = 0

                    for aa in protein_seq:
                        if aa == "-":
                            codon_seq.append("---")
                        else:
                            if nuc_index + 3 <= len(nuc_seq):
                                codon_seq.append(nuc_seq[nuc_index : nuc_index + 3])
                                nuc_index += 3
                            else:
                                break

                    if codon_seq:
                        codon_record = SeqRecord(
                            Seq("".join(codon_seq)),
                            id=protein_record.id,
                            description="",
                        )
                        codon_records.append(codon_record)

                # Write codon alignment
                if codon_records:
                    SeqIO.write(codon_records, output_file, "fasta")

            except Exception as e:
                self.logger.error(f"Failed to create codon alignment for {og_id}: {e}")

        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            futures = [
                executor.submit(create_single_codon_alignment, og_id)
                for og_id in qualifying_orthogroups
            ]

            for future in tqdm(
                as_completed(futures),
                total=len(futures),
                desc="Creating codon alignments",
            ):
                future.result()


@app.command()
def main(
    output_dir: str = typer.Option(
        ..., "-o", "--output-dir", help="Output directory for results"
    ),
    db_query: bool = typer.Option(
        False, "--db-query", help="Query and download genomes from database"
    ),
    fasta_input: Optional[str] = typer.Option(
        None,
        "--fasta-input",
        help="Path to FASTA file or directory containing genomes",
    ),
    fasta_dir: Optional[str] = typer.Option(
        None, "--fasta-dir", help="Directory containing FASTA files"
    ),
    taxonomy: Optional[str] = typer.Option(
        None, "-x", "--taxonomy", help="Taxonomy string for database query"
    ),
    tax_level: str = typer.Option(
        "species", "--tax-level", help="Taxonomic level for grouping"
    ),
    conservation: float = typer.Option(
        0.3, "-c", "--conservation", help="Conservation level for orthogroups"
    ),
    threads: int = typer.Option(4, "-j", "--threads", help="Number of threads"),
    limit: Optional[int] = typer.Option(
        None, "-n", "--limit", help="Limit number of genomes"
    ),
):
    """
    Unified probe design pipeline for bacteria and viruses.

    Two modes of operation:
    1. Database mode (--db-query): Query and download genomes from database
    2. Local mode (--fasta-dir): Use local FASTA files
    """

    # Validate inputs
    if db_query and (fasta_input or fasta_dir):
        typer.echo("Error: Cannot combine --db-query with local FASTA input", err=True)
        raise typer.Exit(1)

    if not db_query and not (fasta_input or fasta_dir):
        typer.echo("Error: Must specify either --db-query or --fasta-input/--fasta-dir", err=True)
        raise typer.Exit(1)

    if db_query and not taxonomy:
        typer.echo("Error: --taxonomy required when using --db-query", err=True)
        raise typer.Exit(1)

    # Initialize pipeline (output directory resolution handled in __init__)
    pipeline = UnifiedProbeDesign(output_dir, threads)

    try:
        if db_query:
            # Database mode
            pipeline.run_database_mode(
                taxonomy=taxonomy,
                tax_level=tax_level,
                limit=limit,
                conservation_level=conservation,
            )
        else:
            # Local mode (file or directory)
            fasta_path = fasta_input if fasta_input else fasta_dir
            if fasta_input is None and fasta_dir is not None:
                typer.echo(
                    "[INFO] --fasta-dir is supported; consider --fasta-input (accepts files or directories).",
                    err=True,
                )
            pipeline.run_local_mode(fasta_input=fasta_path, conservation_level=conservation)

        typer.echo(f"Pipeline completed successfully. Results in: {output_dir}")

    except Exception as e:
        typer.echo(f"Pipeline failed: {e}", err=True)
        raise typer.Exit(1)


if __name__ == "__main__":
    app()


@app.callback(invoke_without_command=True)
def _default(
    ctx: typer.Context,
    output_dir: str = typer.Option(
        None, "-o", "--output-dir", help="Output directory for results"
    ),
    db_query: bool = typer.Option(
        False, "--db-query", help="Query and download genomes from database"
    ),
    fasta_input: Optional[str] = typer.Option(
        None,
        "--fasta-input",
        help="Path to FASTA file or directory containing genomes",
    ),
    fasta_dir: Optional[str] = typer.Option(
        None, "--fasta-dir", help="Directory containing FASTA files"
    ),
    taxonomy: Optional[str] = typer.Option(
        None, "-x", "--taxonomy", help="Taxonomy string for database query"
    ),
    tax_level: str = typer.Option(
        "species", "--tax-level", help="Taxonomic level for grouping"
    ),
    conservation: float = typer.Option(
        0.3, "-c", "--conservation", help="Conservation level for orthogroups"
    ),
    threads: int = typer.Option(4, "-j", "--threads", help="Number of threads"),
    limit: Optional[int] = typer.Option(
        None, "-n", "--limit", help="Limit number of genomes"
    ),
):
    if ctx.invoked_subcommand is None:
        if output_dir is None:
            raise typer.BadParameter("--output-dir is required")
        return main(
            output_dir=output_dir,
            db_query=db_query,
            fasta_input=fasta_input,
            fasta_dir=fasta_dir,
            taxonomy=taxonomy,
            tax_level=tax_level,
            conservation=conservation,
            threads=threads,
            limit=limit,
        )
    return None
