"""Module for querying and downloading genomes from NeLLi genome database."""

import os
import logging
try:
    import duckdb
except ImportError:
    duckdb = None
import pandas as pd
from typing import List, Dict, Optional
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import time

# Configure logging
logger = logging.getLogger(__name__)

# Database location
DEFAULT_DB_PATH = "/clusterfs/jgi/scratch/science/mgs/nelli/databases/nelli-genomes-db/resources/database/gtdb_genomes.duckdb"
DB_PATH = os.environ.get("PPDESIGN_DUCKDB_PATH", DEFAULT_DB_PATH)


class GenomeDatabase:
    """Interface to query and download genomes from NeLLi genome database."""

    def __init__(self, db_path: str = DB_PATH):
        """Initialize database connection.

        Args:
            db_path: Path to DuckDB database file
        """
        if duckdb is None:
            raise ImportError("duckdb is required to query the genome database. Install duckdb or set PPDESIGN_DUCKDB_PATH to override the default.")

        self.db_path = db_path
        if not os.path.exists(self.db_path):
            raise FileNotFoundError(
                f"Genome database not found at {self.db_path}. Set PPDESIGN_DUCKDB_PATH to override the default."
            )

        # Initialize connection
        self.conn = duckdb.connect(self.db_path, read_only=True)
        logger.info(f"Connected to genome database at {self.db_path}")

    def __del__(self):
        """Close database connection."""
        if hasattr(self, "conn"):
            self.conn.close()

    def query_by_taxonomy(
        self,
        taxonomy: str,
        tax_level: Optional[str] = None,
        limit: Optional[int] = None,
        max_genome_size_mb: Optional[float] = None,
    ) -> pd.DataFrame:
        """Query genomes by taxonomy string.

        Args:
            taxonomy: Taxonomy string to search for (e.g., "Faserviricetes", "p__Pseudomonadota")
            tax_level: Taxonomic level to group by (species, genus, family, etc.)
            limit: Maximum number of genomes to return
            max_genome_size_mb: Maximum genome size in MB

        Returns:
            DataFrame with genome information
        """
        # Build query based on taxonomy format
        if taxonomy.startswith(("d__", "p__", "c__", "o__", "f__", "g__", "s__")):
            # GTDB-style taxonomy
            tax_parts = taxonomy.split(";")
            conditions = []
            for part in tax_parts:
                if part:
                    level_prefix = part[:3]
                    level_map = {
                        "d__": "domain",
                        "p__": "phylum",
                        "c__": "class",
                        "o__": "order_name",
                        "f__": "family",
                        "g__": "genus",
                        "s__": "species",
                    }
                    if level_prefix in level_map:
                        col = level_map[level_prefix]
                        value = part[3:]  # Remove prefix
                        conditions.append(f"{col} = '{value}'")

            where_clause = " AND ".join(conditions) if conditions else "1=1"
        else:
            # Search across all taxonomy levels
            where_clause = f"""
            (domain LIKE '%{taxonomy}%' OR 
             phylum LIKE '%{taxonomy}%' OR 
             class LIKE '%{taxonomy}%' OR 
             order_name LIKE '%{taxonomy}%' OR 
             family LIKE '%{taxonomy}%' OR 
             genus LIKE '%{taxonomy}%' OR 
             species LIKE '%{taxonomy}%')
            """

        # Add size filter if specified
        if max_genome_size_mb:
            where_clause += f" AND genome_size_bp <= {max_genome_size_mb * 1_000_000}"

        # Build full query
        query = f"""
        SELECT 
            assembly_accession as accession,
            organism_name,
            domain,
            phylum,
            class,
            order_name as "order",
            family,
            genus,
            species,
            genome_size_bp as genome_size,
            ftp_path,
            ftp_path as file_path
        FROM gtdb_genomes
        WHERE {where_clause}
        """

        # Add grouping if tax_level specified
        if tax_level:
            # Quote the tax_level if it's a reserved word
            tax_level_quoted = f'"{tax_level}"' if tax_level == "order" else tax_level
            query = f"""
            WITH base_query AS ({query}),
            ranked AS (
                SELECT *,
                    ROW_NUMBER() OVER (PARTITION BY {tax_level_quoted} ORDER BY genome_size) as rn
                FROM base_query
            )
            SELECT * FROM ranked WHERE rn = 1
            """

        # Add limit if specified
        if limit:
            query += f" LIMIT {limit}"

        logger.info(f"Executing query for taxonomy: {taxonomy}")
        result_df = self.conn.execute(query).df()
        logger.info(f"Found {len(result_df)} genomes matching criteria")

        return result_df

    def query_bacteria(
        self,
        genus: Optional[str] = None,
        species: Optional[str] = None,
        phylum: Optional[str] = None,
        limit: Optional[int] = None,
    ) -> pd.DataFrame:
        """Convenience method to query bacterial genomes.

        Args:
            genus: Bacterial genus (e.g., "Escherichia")
            species: Bacterial species (e.g., "Escherichia coli")
            phylum: Bacterial phylum (e.g., "Proteobacteria")
            limit: Maximum number of genomes to return

        Returns:
            DataFrame with bacterial genome information
        """
        # Build taxonomy string
        if species:
            # Species query
            taxonomy = f"s__{species}"
        elif genus:
            # Genus query
            taxonomy = f"g__{genus}"
        elif phylum:
            # Phylum query
            taxonomy = f"p__{phylum}"
        else:
            # All bacteria
            taxonomy = "d__Bacteria"

        return self.query_by_taxonomy(
            taxonomy=taxonomy, tax_level="species" if not phylum else None, limit=limit
        )

    def download_genomes(
        self,
        genomes_df: pd.DataFrame,
        output_dir: str,
        download_fna: bool = True,
        download_faa: bool = True,
        threads: int = 5,
        prefix_rules: Optional[Dict] = None,
    ) -> Dict[str, List[str]]:
        """Download genome files for given accessions.

        Args:
            genomes_df: DataFrame with genome information (must have accession and ftp_path columns)
            output_dir: Directory to save downloaded files
            download_fna: Whether to download genome FASTA files
            download_faa: Whether to download protein FASTA files
            threads: Number of parallel downloads
            prefix_rules: Dictionary mapping taxonomy patterns to prefixes

        Returns:
            Dictionary with 'fna' and 'faa' keys containing paths to downloaded files
        """
        # Create output directories
        os.makedirs(output_dir, exist_ok=True)
        fna_dir = os.path.join(output_dir, "fna")
        faa_dir = os.path.join(output_dir, "faa")

        if download_fna:
            os.makedirs(fna_dir, exist_ok=True)
        if download_faa:
            os.makedirs(faa_dir, exist_ok=True)

        # Prepare download tasks
        download_tasks = []
        for _, row in genomes_df.iterrows():
            accession = row["accession"]
            ftp_path = row["ftp_path"]

            # Determine prefix based on taxonomy
            prefix = self._get_prefix(row, prefix_rules)

            # Format accession (replace dots with hyphens)
            accession_formatted = accession.replace(".", "-")

            if download_fna and pd.notna(ftp_path):
                fna_url = f"{ftp_path}/{os.path.basename(ftp_path)}_genomic.fna.gz"
                fna_output = os.path.join(fna_dir, f"{prefix}{accession_formatted}.fna")
                download_tasks.append((fna_url, fna_output, "fna"))

            if download_faa and pd.notna(ftp_path):
                faa_url = f"{ftp_path}/{os.path.basename(ftp_path)}_protein.faa.gz"
                faa_output = os.path.join(faa_dir, f"{prefix}{accession_formatted}.faa")
                download_tasks.append((faa_url, faa_output, "faa"))

        # Download files in parallel
        downloaded_files = {"fna": [], "faa": []}

        with ThreadPoolExecutor(max_workers=threads) as executor:
            future_to_task = {
                executor.submit(self._download_file, url, output, file_type): (
                    url,
                    output,
                    file_type,
                )
                for url, output, file_type in download_tasks
            }

            for future in tqdm(
                as_completed(future_to_task),
                total=len(download_tasks),
                desc="Downloading genomes",
            ):
                url, output, file_type = future_to_task[future]
                try:
                    success = future.result()
                    if success:
                        downloaded_files[file_type].append(output)
                except Exception as e:
                    logger.error(f"Failed to download {url}: {e}")

        # Post-process: generate missing FAA files with prodigal-gv
        if download_faa and downloaded_files["fna"]:
            self._generate_missing_proteins(downloaded_files, fna_dir, faa_dir, threads)

        return downloaded_files

    def _get_prefix(self, row: pd.Series, prefix_rules: Optional[Dict] = None) -> str:
        """Determine prefix for a genome based on taxonomy."""
        if not prefix_rules:
            # Default prefixes
            if row["domain"] == "Bacteria":
                return "BAC__"
            elif row["domain"] == "Archaea":
                return "ARC__"
            elif row["domain"] == "Eukaryota":
                return "EUK__"
            elif row["domain"] == "Viruses":
                return "VIR__"
            else:
                return "UNK__"

        # Apply custom prefix rules
        taxonomy_full = f"{row['domain']};{row['phylum']};{row['class']};{row['order']};{row['family']};{row['genus']};{row['species']}"

        # Check special 'default' key first
        if "default" in prefix_rules:
            default_prefix = prefix_rules["default"]
        else:
            # Use domain-based default
            domain_defaults = {
                "Bacteria": "BAC__",
                "Archaea": "ARC__",
                "Eukaryota": "EUK__",
                "Viruses": "VIR__",
            }
            default_prefix = domain_defaults.get(row["domain"], "UNK__")

        # Check pattern rules
        for pattern, prefix in prefix_rules.items():
            if pattern != "default" and pattern in taxonomy_full:
                return prefix

        return default_prefix

    def _download_file(self, url: str, output_path: str, file_type: str) -> bool:
        """Download a single file with retries."""
        max_retries = 3

        for attempt in range(max_retries):
            try:
                # Download with wget
                cmd = ["wget", "-q", "-O", f"{output_path}.gz", url]
                result = subprocess.run(cmd, capture_output=True, text=True)

                if result.returncode == 0:
                    # Decompress
                    subprocess.run(["gunzip", "-f", f"{output_path}.gz"], check=True)
                    return True
                else:
                    if attempt < max_retries - 1:
                        time.sleep(2**attempt)  # Exponential backoff
                    else:
                        logger.warning(
                            f"Failed to download {url} after {max_retries} attempts"
                        )
                        return False

            except Exception as e:
                logger.error(f"Error downloading {url}: {e}")
                if attempt < max_retries - 1:
                    time.sleep(2**attempt)
                else:
                    return False

        return False

    def _generate_missing_proteins(
        self, downloaded_files: Dict, fna_dir: str, faa_dir: str, threads: int
    ):
        """Generate protein files for genomes missing FAA files using prodigal-gv."""
        # Find FNA files without corresponding FAA files
        fna_basenames = {
            os.path.basename(f).replace(".fna", "") for f in downloaded_files["fna"]
        }
        faa_basenames = {
            os.path.basename(f).replace(".faa", "") for f in downloaded_files["faa"]
        }

        missing_faa = fna_basenames - faa_basenames

        if missing_faa:
            logger.info(
                f"Generating protein files for {len(missing_faa)} genomes using prodigal-gv"
            )

            tasks = []
            for basename in missing_faa:
                fna_path = os.path.join(fna_dir, f"{basename}.fna")
                faa_path = os.path.join(faa_dir, f"{basename}.faa")
                if os.path.exists(fna_path):
                    tasks.append((fna_path, faa_path))

            with ThreadPoolExecutor(max_workers=threads) as executor:
                futures = [
                    executor.submit(self._run_prodigal_gv, fna, faa)
                    for fna, faa in tasks
                ]

                for future in tqdm(
                    as_completed(futures), total=len(tasks), desc="Running prodigal-gv"
                ):
                    try:
                        future.result()
                    except Exception as e:
                        logger.error(f"Prodigal-gv failed: {e}")

    def _run_prodigal_gv(self, fna_path: str, faa_path: str):
        """Run prodigal-gv on a single genome."""
        cmd = ["prodigal-gv", "-i", fna_path, "-a", faa_path, "-p", "meta", "-q"]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"prodigal-gv failed on {fna_path}: {result.stderr}")

        # Format headers to match expected format
        self._format_prodigal_headers(faa_path)

    def _format_prodigal_headers(self, faa_path: str):
        """Format prodigal output headers to match expected format."""
        # Read the file and update headers
        from Bio import SeqIO

        records = []
        genome_id = os.path.basename(faa_path).replace(".faa", "")

        for record in SeqIO.parse(faa_path, "fasta"):
            # Update the ID to include genome ID
            original_id = record.id
            record.id = f"{genome_id}|{original_id}"
            record.description = ""
            records.append(record)

        # Write back
        SeqIO.write(records, faa_path, "fasta")
