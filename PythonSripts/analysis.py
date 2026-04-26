"""
analysis.py

Analyses terminal stop codon usage in bacterial genomes using FASTA genome
files and Bakta JSON annotation files. Only the final stop codon of each CDS
is counted — internal stop codons (which Bakta does not normally annotate)
are intentionally excluded.
"""

import os
import json
from collections import Counter
from Bio import SeqIO


# ── Path configuration ────────────────────────────────────────────────────────

BASE_DIR  = os.path.dirname(os.path.abspath(__file__))
FA_DIR    = os.path.join(BASE_DIR, "DataFiles", "FAFiles")
JSON_DIR  = os.path.join(BASE_DIR, "DataFiles", "JsonFiles")
OUT_FILE        = os.path.join(BASE_DIR, "stop_codon_results.csv")
CODON_OUT_FILE  = os.path.join(BASE_DIR, "codon_usage_results.csv")

VALID_STOP_CODONS = {"TAA", "TAG", "TGA"}

USE_TEST_FILES = True   # True → only process testFile1; False → process all files


# ── Core analysis function ────────────────────────────────────────────────────

def count_stop_codons_json(fasta_path: str, json_path: str) -> dict:
    """
    Count terminal stop codons for all CDS features in one genome.

    Bakta uses 1-based, inclusive coordinates [start, stop].
    The stop codon sits immediately after the annotated CDS on the + strand,
    or immediately before on the – strand.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA genome file.
    json_path : str
        Path to the Bakta JSON annotation file for the same genome.

    Returns
    -------
    dict with keys: TAA, TAG, TGA, total_valid, cds_count, invalid_count,
                    total_cds_length, codon_counts, genome_size, contig_count,
                    gc_content
    """

    # Load all contigs into a dict keyed by sequence ID
    try:
        genome = {rec.id: rec.seq for rec in SeqIO.parse(fasta_path, "fasta")}
    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    if not genome:
        raise ValueError(f"No sequences found in FASTA file: {fasta_path}")

    # Genome-level stats derived directly from the FASTA
    genome_size   = sum(len(seq) for seq in genome.values())
    contig_count  = len(genome)
    gc_count      = sum(seq.upper().count("G") + seq.upper().count("C")
                        for seq in genome.values())
    gc_content    = round(gc_count / genome_size * 100, 4) if genome_size else 0.0

    # Load the Bakta JSON annotation
    try:
        with open(json_path) as fh:
            annotation = json.load(fh)
    except FileNotFoundError:
        raise FileNotFoundError(f"JSON file not found: {json_path}")
    except json.JSONDecodeError as e:
        raise ValueError(f"Malformed JSON in {json_path}: {e}")

    counts = {"TAA": 0, "TAG": 0, "TGA": 0,
              "total_valid": 0, "cds_count": 0, "invalid_count": 0,
              "total_cds_length": 0, "codon_counts": Counter(),
              "genome_size": genome_size, "contig_count": contig_count,
              "gc_content": gc_content}

    for feature in annotation.get("features", []):

        # Only process protein-coding sequences (Bakta uses lowercase "cds")
        if feature.get("type", "").upper() != "CDS":
            continue

        counts["cds_count"] += 1

        # ── Extract coordinates and strand ────────────────────────────────
        # Bakta JSON uses 1-based start, 1-based inclusive stop, and
        # strand as "+" or "-".
        seq_id = feature.get("contig") or feature.get("sequence")
        strand = feature.get("strand")

        try:
            start_1based = int(feature["start"])  # 1-based start of CDS
            end_1based   = int(feature["stop"])   # 1-based inclusive end of CDS
        except (KeyError, TypeError, ValueError):
            counts["invalid_count"] += 1
            continue

        # Validate that the contig exists in the genome FASTA
        if seq_id not in genome:
            counts["invalid_count"] += 1
            continue

        contig = genome[seq_id]
        counts["total_cds_length"] += end_1based - start_1based + 1

        # ── Extract the terminal stop codon ───────────────────────────────
        # Bakta includes the stop codon inside the annotated coordinates:
        #   + strand: stop codon = last 3 bases  → contig[end_1based-3 : end_1based]
        #   - strand: stop codon = first 3 bases (rev-comp) → contig[start_0 : start_0+3].rc
        start_0 = start_1based - 1   # convert to 0-based

        try:
            if strand == "+":
                stop_codon = str(contig[end_1based - 3 : end_1based]).upper()

            elif strand == "-":
                stop_codon = str(contig[start_0 : start_0 + 3]
                                 .reverse_complement()).upper()

            else:
                counts["invalid_count"] += 1
                continue

        except Exception:
            counts["invalid_count"] += 1
            continue

        # ── Count all in-frame codons in this CDS (body + stop codon) ────
        try:
            if strand == "+":
                full_cds = str(contig[start_0 : end_1based]).upper()
            else:
                full_cds = str(contig[start_0 : end_1based]
                               .reverse_complement()).upper()
            for i in range(0, len(full_cds) - 2, 3):
                triplet = full_cds[i : i + 3]
                if len(triplet) == 3:
                    counts["codon_counts"][triplet] += 1
        except Exception:
            pass

        # ── Tally the stop codon ──────────────────────────────────────────
        if stop_codon in VALID_STOP_CODONS:
            counts[stop_codon]  += 1
            counts["total_valid"] += 1
        else:
            # Codon was not TAA, TAG, or TGA — likely a coordinate issue
            counts["invalid_count"] += 1

    return counts


# ── Proportion helper ─────────────────────────────────────────────────────────

def safe_proportion(part: int, total: int) -> str:
    """Return proportion as a string rounded to 4 dp, or 'N/A' if total is 0."""
    if total == 0:
        return "N/A"
    return f"{part / total:.4f}"


# ── File-matching function ────────────────────────────────────────────────────

def match_genome_files(fa_dir: str = FA_DIR,
                       json_dir: str = JSON_DIR
                       ) -> list[tuple[str, str, str]]:
    """
    Pair every 'name.fa' in fa_dir with its 'name.bakta.json' in json_dir.

    Returns a list of (genome_id, fasta_path, json_path) tuples for each
    matched pair.  Unmatched files on either side are reported to stdout.
    """
    fa_map: dict[str, str] = {}
    for fname in sorted(os.listdir(fa_dir)):
        if fname.endswith(".fa"):
            genome_id = fname.split(".")[0]
            fa_map[genome_id] = os.path.join(fa_dir, fname)

    json_map: dict[str, str] = {}
    for fname in sorted(os.listdir(json_dir)):
        if fname.endswith(".json"):
            genome_id = fname.split(".")[0]
            json_map[genome_id] = os.path.join(json_dir, fname)

    matched: list[tuple[str, str, str]] = []
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


# ── Batch scan function ───────────────────────────────────────────────────────

def scan_genome_folder(fa_dir: str = FA_DIR,
                       json_dir: str = JSON_DIR,
                       output_path: str = OUT_FILE,
                       codon_output_path: str = CODON_OUT_FILE) -> None:
    """
    Match FASTA and JSON files by genome ID, run stop-codon analysis on each
    pair, and write a tab-separated results table.

    Genome ID is derived from the filename stem, e.g.
    'SAMEA7491883.fa'  →  genome_id = 'SAMEA7491883'
    'SAMEA7491883.bakta.json'  →  genome_id = 'SAMEA7491883'

    The output file has the following columns:
        Genome, TAA, TAG, TGA, Total_Valid_Stop_Codons,
        CDS_Count, Invalid_Count,
        TAA_Proportion, TAG_Proportion, TGA_Proportion
    """

    pairs = match_genome_files(fa_dir, json_dir)
    if not pairs:
        print("No matched genome/JSON pairs found.")
        return

    header = (
        "Genome,TAA,TAG,TGA,"
        "Total_Valid_Stop_Codons,CDS_Count,Invalid_Count,"
        "TAA_Proportion,TAG_Proportion,TGA_Proportion,Avg_CDS_Length"
    )

    results: list[tuple[str, dict]] = []

    for genome_id, fasta_path, json_path in pairs:
        try:
            c = count_stop_codons_json(fasta_path, json_path)
        except Exception as e:
            print(f"  [ERROR] {genome_id}: {e}")
            continue

        results.append((genome_id, c))
        total = c["total_valid"]
        print(f"  [OK] {genome_id}  —  "
              f"TAA={c['TAA']}, TAG={c['TAG']}, TGA={c['TGA']}, "
              f"Valid={total}, CDS={c['cds_count']}, Invalid={c['invalid_count']}")

    # ── Write stop-codon summary ──────────────────────────────────────────────
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

    # ── Write full codon-usage table ──────────────────────────────────────────
    # Collect all codons seen across every genome for consistent column order
    all_codons = sorted({c for _, res in results for c in res["codon_counts"]})

    with open(codon_output_path, "w") as cout:
        cout.write("Genome," + ",".join(all_codons) + "\n")
        for genome_id, c in results:
            counts_row = ",".join(str(c["codon_counts"].get(cod, 0))
                                  for cod in all_codons)
            cout.write(f"{genome_id},{counts_row}\n")
    print(f"Codon-usage results written to: {codon_output_path}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    scan_genome_folder()
