def has_size_annotations(fasta_path):
    """Check if first sequence header contains ;size=  and ;length= from derep."""
    
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                return ';size=' in line and ';length=' in line
    return False
