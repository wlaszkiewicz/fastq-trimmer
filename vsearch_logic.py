import os
import gzip
import subprocess
from Bio import SeqIO




# def fastq_gz_to_fasta(fastq_gz_path, fasta_path):
#     """Convert a .fastq.gz file to a plain .fasta file."""
#     with gzip.open(fastq_gz_path, "rt") as infile, open(fasta_path, "w") as outfile:
#         for record in SeqIO.parse(infile, "fastq"):
#             outfile.write(f">{record.id}\n{str(record.seq)}\n")


def run_derep(fastq_path, output_path):
    """Run vsearch dereplication on a fasta file."""
   
    cmd = [
        "vsearch",
        "--fastx_uniques", fastq_path,
        "--fastaout", output_path,
        "--fastq_qmax", "60",
        "--sizeout",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"vsearch derep failed:\n{result.stderr}")
    return result.stderr 


def run_fastq_stats(fastq_path, log_path):
    """Run vsearch --fastq_stats and save output to log file."""
    cmd = [
        "vsearch",
        "--fastq_stats", fastq_path,
        "--log", log_path,
        "--fastq_qmax", "60"
    ]
  
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"vsearch fastq_stats failed:\n{result.stderr}")
    return log_path
