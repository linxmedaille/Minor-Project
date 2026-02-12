from Bio import SeqIO
import os

import os
from Bio import SeqIO

# Get folder where script is
base_dir = os.path.dirname(os.path.abspath(__file__))
FAFile = "Mycobacterium_tuberculosis_gca_001196295.6505_4_11_.dna.toplevel.fa"
GffFile = "Mycobacterium_tuberculosis_gca_001196295.6505_4_11.62.gff3"

FATestFile = "test1.fa"
GffTestFile = "test1.gff3"
# Construct full paths
currentFaFile = os.path.join(base_dir, "DataFiles", "FAFiles", FATestFile )
currentGffFile = os.path.join(base_dir, "DataFiles", "GFF3Files", GffTestFile)


# Load genome
genome = {rec.id: rec.seq for rec in SeqIO.parse(currentFaFile, "fasta")}

# Codons that signal stop
stop_codons = {"TAA", "TAG", "TGA"}
stop_count = 0


with open(currentGffFile) as gff:
    TAACount = 0
    TAGCount = 0
    TGACount = 0
    gff.seek(0)
    for line in gff:
        if line.startswith("#"):
            print("skip")
            continue
            
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        
        seq_id, source, feature_type, start, end, score, strand, phase, attributes = parts
        
        if feature_type != "CDS":
            continue
        
        start = int(start) - 1 + int(phase)  # adjust for phase
        end = int(end)
        seq = genome[seq_id][start:end]
        if strand == "-":
            seq = seq.reverse_complement()
        for i in range(0, len(seq) - 2, 3):
            codon = str(seq[i:i+3]).upper()
            if codon in stop_codons:
                match codon:
                    case "TAA":
                        TAACount += 1
                    case "TGA":
                        TGACount += 1
                    case "TAG":
                        TAGCount += 1
                        

print(f"TAA: {TAACount} TGA: {TGACount} TAG: {TAGCount}")