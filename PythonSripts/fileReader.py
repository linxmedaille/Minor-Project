from Bio import SeqIO
import os


# Get project root
base_dir = os.path.dirname(os.path.abspath(__file__))

# Data folders
fa_dir = os.path.join(base_dir, "DataFiles", "FAFiles")
gff_dir = os.path.join(base_dir, "DataFiles", "GFF3Files")
json_path = os.path.join(base_dir, "DataFiles", "JsonFiles")

# Output file
output_file = os.path.join(base_dir, "stop_codon_results.txt")

stop_codons = {"TAA", "TAG", "TGA"}

import json

from Bio import SeqIO
import json

stop_codons = {"TAA", "TAG", "TGA"}

def count_stop_codons_json(fasta_path, json_path):
    genome = {rec.id: rec.seq for rec in SeqIO.parse(fasta_path, "fasta")}

    TAACount = 0
    TAGCount = 0
    TGACount = 0
    codonCount = 0
    cdsCount = 0
    invalidCount = 0

    with open(json_path) as f:
        data = json.load(f)

    for feature in data.get("features", []):
        if feature.get("type") != "CDS":
            continue

        seq_id = feature.get("sequence") or feature.get("contig")
        start = int(feature["start"]) - 1
        end = int(feature.get("stop") or feature.get("end"))
        strand = feature["strand"]

        if seq_id not in genome:
            invalidCount += 1
            continue

        contig = genome[seq_id]
        cdsCount += 1

        if strand == "+":
            stop_codon = str(contig[end:end+3]).upper()
        else:
            stop_codon = str(contig[max(0, start-3):start].reverse_complement()).upper()

        if stop_codon == "TAA":
            TAACount += 1
            codonCount += 1
        elif stop_codon == "TAG":
            TAGCount += 1
            codonCount += 1
        elif stop_codon == "TGA":
            TGACount += 1
            codonCount += 1
        else:
            invalidCount += 1

    return TAACount, TAGCount, TGACount, codonCount, cdsCount, invalidCount



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

            TAA, TAG, TGA, total = count_stop_codons(fasta_path, gff_path)

            out.write(f"{fasta_file}\t{TAA}\t{TAG}\t{TGA}\t{total}\n")

            processed.add(fasta_file)


#scan_files()
#count_stop_codons_json()
#print("Finished")

if __name__ == "__main__":

    fasta_path = os.path.join(fa_dir, "SAMEA7491883.fa")
    json_path = os.path.join(json_path, "SAMEA7491883.bakta.json")

    TAACount, TAGCount, TGACount, codonCount, cdsCount, invalidCount = count_stop_codons_json(fasta_path, json_path)

    print("Results for SAMEA7491883:")
    print(f"TAA: {TAACount}")
    print(f"TAG: {TAGCount}")
    print(f"TGA: {TGACount}")
    print(f"Total codons: {codonCount}")
    