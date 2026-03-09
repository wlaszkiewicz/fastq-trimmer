import os
import glob
import gzip
from Bio.Seq import Seq
from Bio import SeqIO

# references 
ref1     = "TATCGAGAAA"
ref2     = "TTTCAAT"
ref1_rev = str(Seq(ref1).reverse_complement())
ref2_rev = str(Seq(ref2).reverse_complement())


def merge_fastq_gz(barcode_dir, merged_path):
    files = sorted(glob.glob(os.path.join(barcode_dir, "*.fastq.gz")))
    if not files:
        raise FileNotFoundError(f"No .fastq.gz files found in {barcode_dir}")
    with open(merged_path, "wb") as outfile:
        for f in files:
            with open(f, "rb") as infile:
                outfile.write(infile.read())
    return merged_path


def count_records(gz_path):
    count = 0
    with gzip.open(gz_path, "rt") as f:
        for _ in SeqIO.parse(f, "fastq"):
            count += 1
    return count


def _slice_record(record, seq_str, start, end):
    orig_qual = record.letter_annotations["phred_quality"]
    record.letter_annotations = {}
    record.seq = Seq(seq_str[start:end])
    record.letter_annotations["phred_quality"] = orig_qual[start:end]
    return record


def trim_record(record):
    seq = str(record.seq)
    f1  = ref1 in seq
    f2  = ref2 in seq
    r1  = ref1_rev in seq
    r2  = ref2_rev in seq

    if f1 and f2 and not (r1 or r2):
        start = seq.find(ref1)
        end   = seq.find(ref2) + len(ref2)
        return _slice_record(record, seq, start, end), True

    elif r1 and r2 and not (f1 or f2):
        seq_rc = str(record.seq.reverse_complement())
        start  = seq_rc.find(ref1)
        end    = seq_rc.find(ref2) + len(ref2)
        record.letter_annotations["phred_quality"] = record.letter_annotations["phred_quality"][::-1]
        return _slice_record(record, seq_rc, start, end), True

    else:
        return record, False


def process_file(input_path, output_path, incomplete_path, progress_callback=None):
    total = count_records(input_path)
    if total == 0:
        raise ValueError("No records found in input file.")

    records_out        = []
    records_incomplete = []

    with gzip.open(input_path, "rt") as infile:
        for i, record in enumerate(SeqIO.parse(infile, "fastq")):
            trimmed, is_complete = trim_record(record)
            if is_complete:
                records_out.append(trimmed)
            else:
                records_incomplete.append(trimmed)
            if progress_callback:
                progress_callback(int((i + 1) / total * 100))

    with gzip.open(output_path, "wt") as outfile:
        SeqIO.write(records_out, outfile, "fastq")
    with gzip.open(incomplete_path, "wt") as incfile:
        SeqIO.write(records_incomplete, incfile, "fastq")

    return len(records_out), len(records_incomplete)
