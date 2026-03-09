"""
Two layers of testing:
  1. CORRECTNESS  — does the output match the known expected result?
  2. PARITY       — does new fast logic match old BioPython logic exactly?
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ref1     = "TATCGAGAAA"
ref2     = "TTTCAAT"
ref1_rev = str(Seq(ref1).reverse_complement())
ref2_rev = str(Seq(ref2).reverse_complement())


def old_trim_record(record):
    seq = str(record.seq)
    f1  = ref1 in seq
    f2  = ref2 in seq
    r1  = ref1_rev in seq
    r2  = ref2_rev in seq

    if f1 and f2 and not (r1 or r2):
        start     = seq.find(ref1)
        end       = seq.find(ref2) + len(ref2)
        orig_qual = record.letter_annotations["phred_quality"]
        record.letter_annotations = {}
        record.seq = Seq(seq[start:end])
        record.letter_annotations["phred_quality"] = orig_qual[start:end]
        return str(record.seq), record.letter_annotations["phred_quality"], True

    elif r1 and r2 and not (f1 or f2):
        seq_rc    = str(record.seq.reverse_complement())
        start     = seq_rc.find(ref1)
        end       = seq_rc.find(ref2) + len(ref2)
        orig_qual = record.letter_annotations["phred_quality"][::-1]
        record.letter_annotations = {}
        record.seq = Seq(seq_rc[start:end])
        record.letter_annotations["phred_quality"] = orig_qual[start:end]
        return str(record.seq), record.letter_annotations["phred_quality"], True

    else:
        return str(record.seq), record.letter_annotations["phred_quality"], False


_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def _revcomp(seq):
    return seq.translate(_COMP)[::-1]

def new_trim_record(seq, qual_str):
    f1 = ref1 in seq
    f2 = ref2 in seq
    r1 = ref1_rev in seq
    r2 = ref2_rev in seq

    if f1 and f2 and not (r1 or r2):
        start = seq.index(ref1)
        end   = seq.index(ref2) + len(ref2)
        return seq[start:end], qual_str[start:end], True

    elif r1 and r2 and not (f1 or f2):
        seq_rc  = _revcomp(seq)
        qual_rc = qual_str[::-1]
        start   = seq_rc.index(ref1)
        end     = seq_rc.index(ref2) + len(ref2)
        return seq_rc[start:end], qual_rc[start:end], True

    else:
        return seq, qual_str, False


def phred_to_str(phred_list):
    return "".join(chr(q + 33) for q in phred_list)

def str_to_phred(qual_str):
    return [ord(c) - 33 for c in qual_str]

def make_record(seq_str, qual_list=None):
    if qual_list is None:
        qual_list = [30] * len(seq_str)
    return SeqRecord(
        Seq(seq_str),
        id="test",
        description="",
        letter_annotations={"phred_quality": qual_list}
    )

def make_qual_str(length, value=30):
    return phred_to_str([value] * length)


passed = 0
failed = 0

def check(label, condition, detail=""):
    global passed, failed
    if condition:
        print(f"    PASS  {label}")
        passed += 1
    else:
        print(f"    FAIL  {label}{(' — ' + detail) if detail else ''}")
        failed += 1


def run_correctness_tests():
    print("\n" + "="*60)
    print("LAYER 1 — Correctness (new logic vs known expected result)")
    print("="*60)

    payload_fwd = "GGGGGGGGGG"
    payload_rev = "CCCCCCCCCC"

    print("\n[1] Clean forward read")
    seq      = "AAAA" + ref1 + payload_fwd + ref2 + "TTTT"
    qual_str = make_qual_str(len(seq))
    expected_seq  = ref1 + payload_fwd + ref2
    expected_qual = make_qual_str(len(expected_seq))

    result_seq, result_qual, ok = new_trim_record(seq, qual_str)
    check("is complete",            ok)
    check("seq starts with ref1",   result_seq.startswith(ref1))
    check("seq ends with ref2",     result_seq.endswith(ref2))
    check("seq == expected",        result_seq == expected_seq,
          f"got '{result_seq}' expected '{expected_seq}'")
    check("qual length == seq length", len(result_qual) == len(result_seq))
    check("qual == expected",       result_qual == expected_qual)

    print("\n[2] Clean reverse read")
    forward_seq  = ref1 + payload_rev + ref2
    rev_seq      = _revcomp(forward_seq)
    full_seq     = "AAAA" + rev_seq + "TTTT"
    qual_str     = make_qual_str(len(full_seq))

    result_seq, result_qual, ok = new_trim_record(full_seq, qual_str)
    check("is complete",            ok)
    check("seq equals forward",     result_seq == forward_seq,
          f"got '{result_seq}' expected '{forward_seq}'")
    check("seq starts with ref1",   result_seq.startswith(ref1))
    check("seq ends with ref2",     result_seq.endswith(ref2))
    check("qual length == seq length", len(result_qual) == len(result_seq))

    print("\n[3] Quality scores sliced correctly (forward)")
    seq       = "AAAA" + ref1 + "GGG" + ref2 + "TTTT"

    qual_list = list(range(len(seq)))
    qual_str  = phred_to_str(qual_list)

    start = len("AAAA")
    end   = len("AAAA") + len(ref1) + len("GGG") + len(ref2)
    expected_qual = phred_to_str(qual_list[start:end])

    result_seq, result_qual, ok = new_trim_record(seq, qual_str)
    check("is complete",            ok)
    check("qual sliced correctly",  result_qual == expected_qual,
          f"got {str_to_phred(result_qual)} expected {str_to_phred(expected_qual)}")

    print("\n[4] Quality scores sliced correctly (reverse)")
    forward_seq = ref1 + "TTT" + ref2
    rev_seq     = _revcomp(forward_seq)
    full_seq    = "AAAA" + rev_seq + "TTTT"
    qual_list   = list(range(20, 20 + len(full_seq)))
    qual_str    = phred_to_str(qual_list)

    result_seq, result_qual, ok = new_trim_record(full_seq, qual_str)
    check("is complete",            ok)
    check("seq equals forward",     result_seq == forward_seq)
    check("qual length == seq length", len(result_qual) == len(result_seq))

    rev_qual     = qual_str[::-1]
    start        = _revcomp(full_seq).index(ref1)
    end          = _revcomp(full_seq).index(ref2) + len(ref2)
    expected_qual = rev_qual[start:end]
    check("qual reversed+sliced correctly", result_qual == expected_qual)

    print("\n[5] No primers → incomplete")
    seq      = "AAAAAAAAAAAAAAAAAAAAAAAAA"
    qual_str = make_qual_str(len(seq))
    result_seq, result_qual, ok = new_trim_record(seq, qual_str)
    check("is incomplete",          not ok)
    check("seq unchanged",          result_seq == seq)
    check("qual unchanged",         result_qual == qual_str)

    print("\n[6] Only ref1, no ref2 → incomplete")
    seq      = "AAAA" + ref1 + "GGGGGGGG"
    qual_str = make_qual_str(len(seq))
    _, _, ok = new_trim_record(seq, qual_str)
    check("is incomplete",          not ok)

    print("\n[7] Only ref2, no ref1 → incomplete")
    seq      = "GGGGGGGG" + ref2 + "TTTT"
    qual_str = make_qual_str(len(seq))
    _, _, ok = new_trim_record(seq, qual_str)
    check("is incomplete",          not ok)

    print("\n[8] Both forward and reverse primers → incomplete")
    seq      = ref1 + "GGG" + ref2 + "AAA" + ref1_rev + "CCC" + ref2_rev
    qual_str = make_qual_str(len(seq))
    _, _, ok = new_trim_record(seq, qual_str)
    check("is incomplete",          not ok)

    print("\n[9] Reverse complement lookup table")
    check("A→T",   _revcomp("A") == "T")
    check("T→A",   _revcomp("T") == "A")
    check("G→C",   _revcomp("G") == "C")
    check("C→G",   _revcomp("C") == "G")
    check("ATCG→CGAT", _revcomp("ATCG") == "CGAT")
    check("matches BioPython", _revcomp(ref1) == ref1_rev)


def run_parity_tests():
    print("\n" + "="*60)
    print("LAYER 2 — Parity (new fast logic == old BioPython logic)")
    print("="*60)

    cases = [
        ("Clean forward",           "AAAA" + ref1 + "GGGGGGGGGG" + ref2 + "TTTT"),
        ("Clean reverse",           "AAAA" + _revcomp(ref1 + "CCCCCCCCCC" + ref2) + "TTTT"),
        ("No primers",              "AAAAAAAAAAAAAAAAAAAAAAAAA"),
        ("Only ref1",               "AAAA" + ref1 + "GGGGGGGGGG"),
        ("Only ref2",               "GGGGGGGGGG" + ref2 + "TTTT"),
        ("Ambiguous fwd+rev",       ref1 + "GGG" + ref2 + "AAA" + ref1_rev + "CCC" + ref2_rev),
        ("Forward no flanking",     ref1 + "TTTTTTTTTT" + ref2),
        ("Reverse no flanking",     _revcomp(ref1 + "TTTTTTTTTT" + ref2)),
        ("Unique quality scores",   "CCCC" + ref1 + "ATATAT" + ref2 + "GGGG"),
    ]

    for label, seq_str in cases:
        print(f"\n[{label}]")

        qual_list = [20 + (i % 20) for i in range(len(seq_str))]
        qual_str  = phred_to_str(qual_list)

        # old
        record = make_record(seq_str, qual_list[:])
        old_seq, old_qual_list, old_ok = old_trim_record(record)
        old_qual_str = phred_to_str(old_qual_list)

        # new
        new_seq, new_qual_str, new_ok = new_trim_record(seq_str, qual_str)

        check("complete flag matches",   old_ok == new_ok)
        check("trimmed seq matches",     old_seq == new_seq,
              f"\n      old: {old_seq}\n      new: {new_seq}")
        check("trimmed qual matches",    old_qual_str == new_qual_str,
              f"\n      old: {str_to_phred(old_qual_str)}\n      new: {str_to_phred(new_qual_str)}")


def print_summary():
    print(f"\n{'='*60}")
    print(f"TOTAL: {passed} passed, {failed} failed")
    if failed == 0:
        print("All tests passed! New logic is correct and matches old logic.")
    else:
        print("Some tests FAILED — check output above.")
    print("="*60)


if __name__ == "__main__":
    run_correctness_tests()
    run_parity_tests()
    print_summary()