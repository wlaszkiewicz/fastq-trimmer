import re
import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt


def parse_length_distribution(log_path):
    """Extract the Read length distribution table from a vsearch stats log."""
    lengths = []
    counts  = []

    in_table = False
    with open(log_path, "r") as f:
        for line in f:
            if "Read length distribution" in line:
                in_table = True
                continue
            if in_table:
                if line.strip().startswith("L") or line.strip().startswith("-") or line.strip() == "":
                    continue
                if not re.match(r"^\s*(>=\s*\d+|\d+)", line):
                    break
                match = re.match(r"^\s*(?:>=\s*)?(\d+)\s+(\d+)", line)
                if match:
                    lengths.append(int(match.group(1)))
                    counts.append(int(match.group(2)))

    if not lengths:
        raise ValueError("Could not find Read length distribution table in log file.")

    paired = sorted(zip(lengths, counts), key=lambda x: x[0])
    lengths = [p[0] for p in paired]
    counts  = [p[1] for p in paired]

    return lengths, counts


def trim_outliers(lengths, counts, threshold_pct=0.1):
    """
    Drop lengths whose count is below threshold_pct% of total reads.
    Returns the trimmed lists and how many lengths were dropped each side.
    """
    total = sum(counts)
    cutoff = total * (threshold_pct / 100)

    paired = list(zip(lengths, counts))

    # drop sparse entries from the left (short reads)
    while paired and paired[0][1] < cutoff:
        paired.pop(0)

    # drop sparse entries from the right (long reads)
    while paired and paired[-1][1] < cutoff:
        paired.pop()

    if not paired:
        return lengths, counts  #

    return [p[0] for p in paired], [p[1] for p in paired]


def plot_length_distribution(log_path, threshold_pct=0.1):
    """Parse log, plot bar chart, save to PNG, and return the figure."""
    lengths, counts = parse_length_distribution(log_path)
    total_lengths   = len(lengths)

 
    lengths_trimmed, counts_trimmed = trim_outliers(lengths, counts, threshold_pct)
    n_dropped = total_lengths - len(lengths_trimmed)

    fig, ax = plt.subplots(figsize=(12, 5))

    ax.bar(lengths_trimmed, counts_trimmed, width=1.0, color="steelblue", edgecolor="none")

    ax.set_xlabel("Read length (bp)")
    ax.set_ylabel("Number of reads")
    ax.set_title("Read Length Distribution")
    ax.margins(x=0.01)

    length_range = lengths_trimmed[-1] - lengths_trimmed[0]
    # note about trimmed outliers
    step = max(1, len(lengths_trimmed) // 20)
    note = f"x-axis shows every {step}bp step (too many ticks to display all)" if length_range > 60 else ""
    if n_dropped > 0:
        ax.annotate(
            f"{n_dropped} sparse length(s) hidden (< {threshold_pct}% of reads)\n{note}",
            xy=(0.01, 0.97), xycoords="axes fraction",
            fontsize=8, color="gray", va="top"
        )

   
    if length_range <= 60:
        ax.set_xticks(lengths_trimmed)
        ax.set_xticklabels([str(l) for l in lengths_trimmed], rotation=45, ha="right", fontsize=7)
    else:
        tick_positions = lengths_trimmed[::step]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels([str(l) for l in tick_positions], rotation=45, ha="right")

    plt.tight_layout()

    return fig

