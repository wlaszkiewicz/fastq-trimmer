import os
import glob
import gzip
import multiprocessing as mp
from Bio.Seq import Seq

def merge_fastq_gz(barcode_dir, merged_path):
    files = sorted(glob.glob(os.path.join(barcode_dir, "*.fastq.gz")))
    if not files:
        raise FileNotFoundError(f"No .fastq.gz files found in {barcode_dir}")
    with open(merged_path, "wb") as outfile:
        for f in files:
            with open(f, "rb") as infile:
                outfile.write(infile.read())
    return merged_path


def _iter_fastq_raw(fileobj):
    """Yield (header, seq, qual) tuples from an open text fileobj."""
    while True:
        header = fileobj.readline()
        if not header:
            break
        seq = fileobj.readline().rstrip("\n")
        fileobj.readline()  # + line, skip
        qual = fileobj.readline().rstrip("\n")
        if header and seq and qual:
            yield header.rstrip("\n"), seq, qual


def estimate_record_count(gz_path):
    """
    Estimate total records from a small sample + file size.
    Avoids reading the whole file just for a count.
    """
    file_size = os.path.getsize(gz_path)
    sample_records = 0
    sample_bytes = 0

    with gzip.open(gz_path, "rt") as f:
        for header, seq, qual in _iter_fastq_raw(f):
            sample_bytes += len(header) + len(seq) + len(qual) + 3
            sample_records += 1
            if sample_records >= 500:
                break

    if sample_records == 0:
        return 1

    avg_bytes_per_record = sample_bytes / sample_records
    estimated_uncompressed = file_size * 3   # gzip ~25-30% compression ratio
    return max(1, int(estimated_uncompressed / avg_bytes_per_record))


# lookup table 
_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def _revcomp(seq):
    return seq.translate(_COMP)[::-1]


def _trim_record_raw(header, seq, qual, ref1, ref2, ref1_rev, ref2_rev):
    """Returns (header, trimmed_seq, trimmed_qual, is_complete)."""
    f1 = ref1 in seq
    f2 = ref2 in seq
    r1 = ref1_rev in seq
    r2 = ref2_rev in seq

    if f1 and f2 and not (r1 or r2):
        start = seq.index(ref1)
        end = seq.index(ref2) + len(ref2)
        return header, seq[start:end], qual[start:end], True

    elif r1 and r2 and not (f1 or f2):
        seq_rc = _revcomp(seq)
        qual_rc = qual[::-1]
        start = seq_rc.index(ref1)
        end = seq_rc.index(ref2) + len(ref2)
        return header, seq_rc[start:end], qual_rc[start:end], True

    else:
        return header, seq, qual, False


def _process_chunk(args):
    """
    chunk: list of (header, seq, qual)
    Returns (complete_lines, incomplete_lines) as lists of fastq record strings.
    """
    chunk, ref1, ref2, ref1_rev, ref2_rev = args
    complete = []
    incomplete = []
    for header, seq, qual in chunk:
        h, s, q, ok = _trim_record_raw(header, seq, qual, ref1, ref2, ref1_rev, ref2_rev)
        record_str = f"{h}\n{s}\n+\n{q}\n"
        if ok:
            complete.append(record_str)
        else:
            incomplete.append(record_str)
    return complete, incomplete


def process_file(input_path, output_path, incomplete_path, ref1, ref2, progress_callback=None):
    ref1_rev = str(Seq(ref1).reverse_complement())
    ref2_rev = str(Seq(ref2).reverse_complement())
    n_cores = max(1, mp.cpu_count() - 1)  # leave one core free for the UI
    chunk_size = 5000                     # records per worker chunk

    estimated_total = estimate_record_count(input_path)
    complete_count = 0
    incomplete_count = 0
    processed = 0

    with gzip.open(input_path, "rt") as infile, \
         gzip.open(output_path, "wt") as outfile, \
         gzip.open(incomplete_path, "wt") as incfile, \
         mp.Pool(processes=n_cores) as pool:

        chunk = []
        jobs = []

        for header, seq, qual in _iter_fastq_raw(infile):
            chunk.append((header, seq, qual))

            if len(chunk) >= chunk_size:
                jobs.append(pool.apply_async(_process_chunk, ((chunk, ref1, ref2, ref1_rev, ref2_rev),)))
                chunk = []

            # drain finished jobs to keep memory under control
            while len(jobs) >= n_cores * 2:
                complete, incomplete = jobs.pop(0).get()
                outfile.writelines(complete)
                incfile.writelines(incomplete)
                complete_count += len(complete)
                incomplete_count += len(incomplete)
                processed += len(complete) + len(incomplete)
                if progress_callback:
                    progress_callback(min(99, int(processed / estimated_total * 100)))

        # submit any remaining records
        if chunk:
            jobs.append(pool.apply_async(_process_chunk, ((chunk, ref1, ref2, ref1_rev, ref2_rev),)))

        # drain all remaining jobs
        for job in jobs:
            complete, incomplete = job.get()
            outfile.writelines(complete)
            incfile.writelines(incomplete)
            complete_count += len(complete)
            incomplete_count += len(incomplete)
            processed += len(complete) + len(incomplete)
            if progress_callback:
                progress_callback(min(99, int(processed / estimated_total * 100)))

    if progress_callback:
        progress_callback(100)

    return complete_count, incomplete_count