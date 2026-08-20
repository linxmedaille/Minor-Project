import os
import json

from Bio import SeqIO


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FA_DIR = os.path.join(BASE_DIR, "DataFiles", "FAFiles")
JSON_DIR = os.path.join(BASE_DIR, "DataFiles", "JsonFiles")
OUT_FILE = os.path.join(BASE_DIR, "stop_codon_results.csv")
CODON_OUT_FILE = os.path.join(BASE_DIR, "codon_usage_results.csv")

VALID_STOP_CODONS = ["TAA", "TAG", "TGA"]

USE_TEST_FILES = True


def read_fasta(fasta_path):
    genome = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        genome[record.id] = record.seq
    return genome


def count_gc(genome):
    gc_count = 0
    for sequence in genome.values():
        text = str(sequence).upper()
        gc_count += text.count("G")
        gc_count += text.count("C")
    return gc_count


def count_stop_codons_json(fasta_path, json_path):
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    genome = read_fasta(fasta_path)
    if len(genome) == 0:
        raise ValueError(f"No sequences found in FASTA file: {fasta_path}")

    genome_size = 0
    for sequence in genome.values():
        genome_size += len(sequence)

    gc_count = count_gc(genome)
    if genome_size > 0:
        gc_content = round(gc_count / genome_size * 100, 4)
    else:
        gc_content = 0.0

    if not os.path.exists(json_path):
        raise FileNotFoundError(f"JSON file not found: {json_path}")

    with open(json_path) as json_file:
        annotation = json.load(json_file)

    counts = {
        "TAA": 0,
        "TAG": 0,
        "TGA": 0,
        "total_valid": 0,
        "cds_count": 0,
        "invalid_count": 0,
        "total_cds_length": 0,
        "codon_counts": {},
        "genome_size": genome_size,
        "contig_count": len(genome),
        "gc_content": gc_content,
    }

    features = annotation.get("features", [])
    for feature in features:
        if feature.get("type", "").upper() != "CDS":
            continue

        counts["cds_count"] += 1

        contig_name = feature.get("contig")
        if contig_name is None:
            contig_name = feature.get("sequence")

        strand = feature.get("strand")

        if contig_name not in genome:
            counts["invalid_count"] += 1
            continue
        if "start" not in feature or "stop" not in feature:
            counts["invalid_count"] += 1
            continue
        if strand != "+" and strand != "-":
            counts["invalid_count"] += 1
            continue

        contig = genome[contig_name]

        start = int(feature["start"])
        end = int(feature["stop"])
        start_index = start - 1

        counts["total_cds_length"] += end - start + 1

        if strand == "+":
            gene = str(contig[start_index:end]).upper()
        else:
            gene = str(contig[start_index:end].reverse_complement()).upper()

        stop_codon = gene[-3:]

        for offset in range(0, len(gene) - 2, 3):
            codon = gene[offset:offset + 3]
            if codon not in counts["codon_counts"]:
                counts["codon_counts"][codon] = 0
            counts["codon_counts"][codon] += 1

        if stop_codon in VALID_STOP_CODONS:
            counts[stop_codon] += 1
            counts["total_valid"] += 1
        else:
            counts["invalid_count"] += 1

    return counts


def safe_proportion(part, total):
    if total == 0:
        return "N/A"
    return f"{part / total:.4f}"


def find_files(folder, ending):
    file_map = {}
    for filename in sorted(os.listdir(folder)):
        if filename.endswith(ending):
            genome_id = filename.split(".")[0]
            file_map[genome_id] = os.path.join(folder, filename)
    return file_map


def match_genome_files(fa_dir=FA_DIR, json_dir=JSON_DIR):
    fasta_files = find_files(fa_dir, ".fa")
    json_files = find_files(json_dir, ".json")

    matched = []
    for genome_id in sorted(fasta_files):
        if USE_TEST_FILES and genome_id != "testFile1":
            continue
        if genome_id in json_files:
            matched.append((genome_id, fasta_files[genome_id], json_files[genome_id]))
        else:
            print(f"  [SKIP] No JSON for {genome_id}.fa")

    for genome_id in sorted(json_files):
        if USE_TEST_FILES and genome_id != "testFile1":
            continue
        if genome_id not in fasta_files:
            print(f"  [SKIP] No FASTA for {genome_id}.bakta.json")

    return matched


def scan_genome_folder(fa_dir=FA_DIR, json_dir=JSON_DIR,
                       output_path=OUT_FILE, codon_output_path=CODON_OUT_FILE):

    pairs = match_genome_files(fa_dir, json_dir)
    if len(pairs) == 0:
        print("No matched genome/JSON pairs found.")
        return

    results = []
    for genome_id, fasta_path, json_path in pairs:
        try:
            counts = count_stop_codons_json(fasta_path, json_path)
        except Exception as error:
            print(f"  [ERROR] {genome_id}: {error}")
            continue

        results.append((genome_id, counts))
        print(f"  [OK] {genome_id}  -  TAA={counts['TAA']}, TAG={counts['TAG']}, "
              f"TGA={counts['TGA']}, Valid={counts['total_valid']}, "
              f"CDS={counts['cds_count']}, Invalid={counts['invalid_count']}")

    header = ("Genome,TAA,TAG,TGA,"
              "Total_Valid_Stop_Codons,CDS_Count,Invalid_Count,"
              "TAA_Proportion,TAG_Proportion,TGA_Proportion,Avg_CDS_Length")

    with open(output_path, "w") as out_file:
        out_file.write(header + "\n")
        for genome_id, counts in results:
            total = counts["total_valid"]
            taa_share = safe_proportion(counts["TAA"], total)
            tag_share = safe_proportion(counts["TAG"], total)
            tga_share = safe_proportion(counts["TGA"], total)
            avg_length = safe_proportion(counts["total_cds_length"], counts["cds_count"])

            row = (f"{genome_id},"
                   f"{counts['TAA']},{counts['TAG']},{counts['TGA']},"
                   f"{total},{counts['cds_count']},{counts['invalid_count']},"
                   f"{taa_share},{tag_share},{tga_share},{avg_length}")
            out_file.write(row + "\n")

    print(f"\nStop-codon results written to: {output_path}")

    all_codons = []
    for genome_id, counts in results:
        for codon in counts["codon_counts"]:
            if codon not in all_codons:
                all_codons.append(codon)
    all_codons.sort()

    with open(codon_output_path, "w") as codon_file:
        codon_file.write("Genome," + ",".join(all_codons) + "\n")
        for genome_id, counts in results:
            numbers = []
            for codon in all_codons:
                numbers.append(str(counts["codon_counts"].get(codon, 0)))
            codon_file.write(genome_id + "," + ",".join(numbers) + "\n")

    print(f"Codon-usage results written to: {codon_output_path}")


if __name__ == "__main__":
    scan_genome_folder()
