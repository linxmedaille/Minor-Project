import os
import json
from collections import Counter
from Bio import SeqIO


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FA_DIR = os.path.join(BASE_DIR, "DataFiles", "FAFiles")
JSON_DIR = os.path.join(BASE_DIR, "DataFiles", "JsonFiles")
OUT_FILE = os.path.join(BASE_DIR, "stop_codon_results.csv")
CODON_OUT_FILE = os.path.join(BASE_DIR, "codon_usage_results.csv")

VALID_STOP_CODONS = {"TAA", "TAG", "TGA"}

USE_TEST_FILES = True


def count_stop_codons_json(fasta_path, json_path):
    try:
        genome = {rec.id: rec.seq for rec in SeqIO.parse(fasta_path, "fasta")}
    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    if not genome:
        raise ValueError(f"No sequences found in FASTA file: {fasta_path}")

    genome_size = sum(len(seq) for seq in genome.values())
    contig_count = len(genome)
    gc_count = sum(seq.upper().count("G") + seq.upper().count("C") for seq in genome.values())
    gc_content = round(gc_count / genome_size * 100, 4) if genome_size else 0.0

    try:
        with open(json_path) as fh:
            annotation = json.load(fh)
    except FileNotFoundError:
        raise FileNotFoundError(f"JSON file not found: {json_path}")
    except json.JSONDecodeError as e:
        raise ValueError(f"Malformed JSON in {json_path}: {e}")

    counts = {
        "TAA": 0, "TAG": 0, "TGA": 0,
        "total_valid": 0, "cds_count": 0, "invalid_count": 0,
        "total_cds_length": 0, "codon_counts": Counter(),
        "genome_size": genome_size, "contig_count": contig_count,
        "gc_content": gc_content,
    }

    for feature in annotation.get("features", []):
        if feature.get("type", "").upper() != "CDS":
            continue

        counts["cds_count"] += 1

        seq_id = feature.get("contig") or feature.get("sequence")
        strand = feature.get("strand")

        try:
            start_1based = int(feature["start"])
            end_1based   = int(feature["stop"])
        except (KeyError, TypeError, ValueError):
            counts["invalid_count"] += 1
            continue

        if seq_id not in genome:
            counts["invalid_count"] += 1
            continue

        contig = genome[seq_id]
        start_0 = start_1based - 1
        counts["total_cds_length"] += end_1based - start_1based + 1

        # Bakta includes the stop codon inside the annotated coordinates:
        # + strand: last 3 bases;  - strand: first 3 bases (rev-comp)
        try:
            if strand == "+":
                stop_codon = str(contig[end_1based - 3 : end_1based]).upper()
            elif strand == "-":
                stop_codon = str(contig[start_0 : start_0 + 3].reverse_complement()).upper()
            else:
                counts["invalid_count"] += 1
                continue
        except Exception:
            counts["invalid_count"] += 1
            continue

        try:
            if strand == "+":
                full_cds = str(contig[start_0 : end_1based]).upper()
            else:
                full_cds = str(contig[start_0 : end_1based].reverse_complement()).upper()
            for i in range(0, len(full_cds) - 2, 3):
                triplet = full_cds[i : i + 3]
                if len(triplet) == 3:
                    counts["codon_counts"][triplet] += 1
        except Exception:
            pass

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


def match_genome_files(fa_dir=FA_DIR, json_dir=JSON_DIR):
    fa_map = {}
    for fname in sorted(os.listdir(fa_dir)):
        if fname.endswith(".fa"):
            genome_id = fname.split(".")[0]
            fa_map[genome_id] = os.path.join(fa_dir, fname)

    json_map = {}
    for fname in sorted(os.listdir(json_dir)):
        if fname.endswith(".json"):
            genome_id = fname.split(".")[0]
            json_map[genome_id] = os.path.join(json_dir, fname)

    matched = []
    for genome_id, fasta_path in sorted(fa_map.items()):
        if USE_TEST_FILES and genome_id != "testFile1":
            continue
        if genome_id in json_map:
            matched.append((genome_id, fasta_path, json_map[genome_id]))
        else:
            print(f"  [SKIP] No JSON for {genome_id}.fa")

    for genome_id in sorted(json_map):
        if USE_TEST_FILES and genome_id != "testFile1":
            continue
        if genome_id not in fa_map:
            print(f"  [SKIP] No FASTA for {genome_id}.bakta.json")

    return matched


def scan_genome_folder(fa_dir=FA_DIR, json_dir=JSON_DIR, output_path=OUT_FILE, codon_output_path=CODON_OUT_FILE):
    pairs = match_genome_files(fa_dir, json_dir)
    if not pairs:
        print("No matched genome/JSON pairs found.")
        return

    header = (
        "Genome,TAA,TAG,TGA,"
        "Total_Valid_Stop_Codons,CDS_Count,Invalid_Count,"
        "TAA_Proportion,TAG_Proportion,TGA_Proportion,Avg_CDS_Length"
    )

    results = []

    for genome_id, fasta_path, json_path in pairs:
        try:
            c = count_stop_codons_json(fasta_path, json_path)
        except Exception as e:
            print(f"  [ERROR] {genome_id}: {e}")
            continue

        results.append((genome_id, c))
        total = c["total_valid"]
        print(f"  [OK] {genome_id}  —  TAA={c['TAA']}, TAG={c['TAG']}, TGA={c['TGA']}, "
              f"Valid={total}, CDS={c['cds_count']}, Invalid={c['invalid_count']}")

    with open(output_path, "w") as out:
        out.write(header + "\n")
        for genome_id, c in results:
            total = c["total_valid"]
            avg_len = safe_proportion(c["total_cds_length"], c["cds_count"])
            row = (
                f"{genome_id},"
                f"{c['TAA']},{c['TAG']},{c['TGA']},"
                f"{total},{c['cds_count']},{c['invalid_count']},"
                f"{safe_proportion(c['TAA'], total)},"
                f"{safe_proportion(c['TAG'], total)},"
                f"{safe_proportion(c['TGA'], total)},"
                f"{avg_len}"
            )
            out.write(row + "\n")
    print(f"\nStop-codon results written to: {output_path}")

    all_codons = sorted({c for _, res in results for c in res["codon_counts"]})
    with open(codon_output_path, "w") as cout:
        cout.write("Genome," + ",".join(all_codons) + "\n")
        for genome_id, c in results:
            counts_row = ",".join(str(c["codon_counts"].get(cod, 0)) for cod in all_codons)
            cout.write(f"{genome_id},{counts_row}\n")
    print(f"Codon-usage results written to: {codon_output_path}")


if __name__ == "__main__":
    scan_genome_folder()
