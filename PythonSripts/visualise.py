import os

import matplotlib.pyplot as plt
import pandas as pd

from analysis import OUT_FILE


FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Figures")
#drawing 24 cut points across the range
GC_BIN_COUNT = 25 

GENOME_SIZE_BIN_COUNT = 25

STOP_CODON_COLORS = {"TAA": "#2ca02c", "TAG": "#1f3b57", "TGA": "#d62728"}


def load_results(path=OUT_FILE):
    table = pd.read_csv(path)
    return table[table["Total_Valid_Stop_Codons"] > 0].copy()


def average(numbers):
    return sum(numbers) / len(numbers)


def standard_error(numbers):
    mean = average(numbers)

    total = 0.0
    for value in numbers:
        total += (value - mean) ** 2
    standard_deviation = (total / len(numbers)) ** 0.5

    return standard_deviation / len(numbers) ** 0.5


def binned_trend(x_values, y_values, bin_count):
    x_values = list(x_values)
    y_values = list(y_values)

    smallest = min(x_values)
    largest = max(x_values)
    bin_width = (largest - smallest) / bin_count
    if bin_width == 0:
        bin_width = 1

    bins = []
    for i in range(bin_count):
        bins.append([])

    for i in range(len(x_values)):
        bin_index = int((x_values[i] - smallest) / bin_width)
        if bin_index >= bin_count:
            bin_index = bin_count - 1
        bins[bin_index].append(y_values[i])

    centers = []
    means = []
    errors = []
    for i in range(bin_count):
        values = bins[i]
        if len(values) == 0:
            continue
        centers.append(smallest + bin_width * (i + 0.5))
        means.append(average(values))
        errors.append(standard_error(values))

    return centers, means, errors


def plot_gc_vs_stop_codon(table, out_path):
    figure, (left_axis, right_axis) = plt.subplots(1, 2, figsize=(12, 5))

    for codon in STOP_CODON_COLORS:
        color = STOP_CODON_COLORS[codon]
        centers, means, errors = binned_trend(table["GC_Content"],
                                              table[codon + "_Proportion"],
                                              GC_BIN_COUNT)

        lower = []
        upper = []
        for i in range(len(means)):
            lower.append(means[i] - errors[i])
            upper.append(means[i] + errors[i])

        left_axis.plot(centers, means, label=codon, color=color)
        left_axis.fill_between(centers, lower, upper, color=color, alpha=0.2)

    left_axis.set_title(f"All genomes (n = {len(table):,})")
    left_axis.set_xlabel("Genomic GC content (%)")
    left_axis.set_ylabel("Stop codon proportion")
    left_axis.set_ylim(0, 1)
    left_axis.legend()

    columns = ["GC_Content", "TAA_Proportion", "TAG_Proportion", "TGA_Proportion"]
    species_means = table.groupby("Species")[columns].mean()
    species_means["Genome_Count"] = table.groupby("Species").size()

    for codon in STOP_CODON_COLORS:
        right_axis.scatter(species_means["GC_Content"],
                           species_means[codon + "_Proportion"],
                           s=species_means["Genome_Count"] ** 0.5 * 3,
                           color=STOP_CODON_COLORS[codon],
                           alpha=0.7, edgecolor="none", label=codon)

    right_axis.set_title(f"Species means (n = {len(species_means)})")
    right_axis.set_xlabel("Genomic GC content (%)")
    right_axis.set_ylim(0, 1)
    right_axis.legend()

    figure.suptitle("Stop codon proportion vs genomic GC content")
    figure.tight_layout()
    figure.savefig(out_path, dpi=150)
    plt.close(figure)


def plot_correlation_heatmap(table, out_path):
    columns = ["TAA_Proportion", "TAG_Proportion", "TGA_Proportion", "GC_Content",
               "Genome_Size_bp", "CDS_Count", "Avg_CDS_Length", "Contig_Count"]
    labels = ["TAA", "TAG", "TGA", "GC%",
              "Genome size", "CDS count", "Avg CDS len", "Contig count"]

    correlation = table[columns].corr()
    correlation.columns = labels
    correlation.index = labels

    figure, axis = plt.subplots(figsize=(8, 7))
    image = axis.imshow(correlation, cmap="RdBu_r", vmin=-1, vmax=1)
    figure.colorbar(image, ax=axis, label="Pearson r")

    axis.set_xticks(range(len(labels)), labels, rotation=45, ha="right")
    axis.set_yticks(range(len(labels)), labels)

    for row in range(len(labels)):
        for column in range(len(labels)):
            value = correlation.iloc[row, column]
            axis.text(column, row, f"{value:.2f}", ha="center", va="center")

    axis.set_title(f"Correlation matrix - genomic variables (n = {len(table):,} genomes)")
    figure.tight_layout()
    figure.savefig(out_path, dpi=150)
    plt.close(figure)


def plot_genome_size_vs_stop_codon(table, out_path):
    genome_size_mbp = table["Genome_Size_bp"] / 1000000

    codons = list(STOP_CODON_COLORS)
    figure, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)

    for i in range(len(codons)):
        codon = codons[i]
        axis = axes[i]
        proportion = table[codon + "_Proportion"]

        r = genome_size_mbp.corr(proportion)

        axis.scatter(genome_size_mbp, proportion, s=4, alpha=0.15,
                     color=STOP_CODON_COLORS[codon])
        centers, means, errors = binned_trend(genome_size_mbp, proportion,
                                              GENOME_SIZE_BIN_COUNT)
        axis.plot(centers, means, color="black")

        axis.set_title(f"{codon}  (r = {r:.3f})")
        axis.set_xlabel("Genome size (Mbp)")

    axes[0].set_ylabel("Proportion")
    axes[0].set_ylim(0, 1)

    figure.suptitle("Genome size vs stop codon proportion")
    figure.tight_layout()
    figure.savefig(out_path, dpi=150)
    plt.close(figure)


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    table = load_results()

    plot_gc_vs_stop_codon(table, os.path.join(FIGURES_DIR, "figure2_gc_vs_stop_codon.png"))
    plot_correlation_heatmap(table, os.path.join(FIGURES_DIR, "figure3_correlation_heatmap.png"))
    plot_genome_size_vs_stop_codon(table, os.path.join(FIGURES_DIR, "figure4_genome_size_vs_stop_codon.png"))

    print(f"Figures written to: {FIGURES_DIR}")


if __name__ == "__main__":
    main()
