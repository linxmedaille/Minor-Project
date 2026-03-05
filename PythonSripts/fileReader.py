from Bio import SeqIO
import os


# Get project root
base_dir = os.path.dirname(os.path.abspath(__file__))

# Data folders
fa_dir = os.path.join(base_dir, "DataFiles", "FAFiles")
gff_dir = os.path.join(base_dir, "DataFiles", "GFF3Files")

# Output file
output_file = os.path.join(base_dir, "stop_codon_results.txt")

stop_codons = {"TAA", "TAG", "TGA"}


def count_stop_codons(fasta_path, gff_path):

    genome = {rec.id: rec.seq for rec in SeqIO.parse(fasta_path, "fasta")}

    TAACount = 0
    TAGCount = 0
    TGACount = 0
    codonCount = 0

    with open(gff_path) as gff:

        for line in gff:

            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")

            if len(parts) < 9:
                continue

            seq_id, source, feature_type, start, end, score, strand, phase, attributes = parts

            if feature_type != "CDS":
                continue

            start = int(start) - 1 + int(phase)
            end = int(end)

            seq = genome[seq_id][start:end]

            if strand == "-":
                seq = seq.reverse_complement()

            for i in range(0, len(seq) - 2, 3):

                codon = str(seq[i:i+3]).upper()

                if codon in stop_codons:

                    if codon == "TAA":
                        TAACount += 1
                    elif codon == "TAG":
                        TAGCount += 1
                    elif codon == "TGA":
                        TGACount += 1

                    codonCount += 1

    return TAACount, TAGCount, TGACount, codonCount


def scan_files():

    fasta_files = sorted([f for f in os.listdir(fa_dir) if f.endswith(".fa")])
    gff_files = sorted([f for f in os.listdir(gff_dir) if f.endswith(".gff3")])

    processed = set()

    with open(output_file, "w") as out:

        out.write("Genome\tTAA\tTAG\tTGA\tTotal\n")

        for fasta_file, gff_file in zip(fasta_files, gff_files):

            if fasta_file in processed:
                continue

            fasta_path = os.path.join(fa_dir, fasta_file)
            gff_path = os.path.join(gff_dir, gff_file)

            print("Processing:", fasta_file, "with", gff_file)

            TAA, TAG, TGA, total = count_stop_codons(fasta_path, gff_path)

            out.write(f"{fasta_file}\t{TAA}\t{TAG}\t{TGA}\t{total}\n")

            processed.add(fasta_file)


scan_files()

print("Finished. Results written to:", output_file)