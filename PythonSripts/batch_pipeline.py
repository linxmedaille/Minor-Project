"""
batch_pipeline.py

Downloads genome assemblies and Bakta annotations from the AllTheBacteria (ATB)
dataset one batch at a time, runs stop-codon analysis on each sample, appends
results to the shared CSV, and deletes all temporary files before moving to the
next batch — keeping disk usage to roughly one batch archive pair at a time.

Progress is written to pipeline_progress.json after every sample, so the run
can be interrupted and restarted without reprocessing completed samples.
"""

import csv
import json
import os
import tarfile
from collections import defaultdict

import requests

from analysis import count_stop_codons_json, safe_proportion, OUT_FILE

# ── Path configuration ────────────────────────────────────────────────────────

BASE_DIR        = os.path.dirname(os.path.abspath(__file__))
FA_LIST_PATH    = os.path.join(BASE_DIR, "DataFiles", "Lists", "file_list.r0.2.v2.tsv")
BAKTA_LIST_PATH = os.path.join(BASE_DIR, "DataFiles", "Lists", "atb.bakta.r0.2.status.tsv")
PROGRESS_FILE   = os.path.join(BASE_DIR, "pipeline_progress.json")
TEMP_DIR        = os.path.join(BASE_DIR, "temp_downloads")

OSF_API        = "https://api.osf.io/v2"
BAKTA_OSF_NODE = "tyw72"    # OSF project node containing the bakta batch archives


# ── OSF URL discovery ─────────────────────────────────────────────────────────

def _osf_list_files(storage_url: str) -> dict[str, str]:
    """
    Recursively page through an OSF storage URL and return all
    .tar.xz files as filename → download URL, including files in subfolders.
    """
    found: dict[str, str] = {}
    url: str | None = storage_url

    while url:
        r = requests.get(url, params={"page[size]": 100}, timeout=30)
        r.raise_for_status()
        payload = r.json()

        for item in payload["data"]:
            attrs = item["attributes"]
            name  = attrs["name"]
            kind  = attrs["kind"]
            if kind == "file" and name.endswith(".tar.xz"):
                found[name] = item["links"]["download"]
            elif kind == "folder":
                sub_url = item["relationships"]["files"]["links"]["related"]["href"]
                found.update(_osf_list_files(sub_url))

        url = payload["links"].get("next")

    return found


def _fetch_bakta_urls_from_osf() -> dict[str, str]:
    """
    Discover all bakta batch archive download URLs from the ATB OSF component.
    Returns a mapping of archive filename → direct download URL.
    """
    storage_url = f"{OSF_API}/nodes/{BAKTA_OSF_NODE}/files/osfstorage/"
    all_files   = _osf_list_files(storage_url)
    return {name: url for name, url in all_files.items() if "bakta" in name}


# ── Progress tracking ─────────────────────────────────────────────────────────

def _load_progress() -> set:
    if os.path.exists(PROGRESS_FILE):
        with open(PROGRESS_FILE) as f:
            content = f.read().strip()
            if content:
                return set(json.loads(content))
    return set()


def _save_progress(done: set) -> None:
    with open(PROGRESS_FILE, "w") as f:
        json.dump(sorted(done), f, indent=2)


# ── Download helper ───────────────────────────────────────────────────────────

def _download(url: str, dest: str) -> None:
    name = os.path.basename(dest)
    print(f"  Downloading {name} ...", end="", flush=True)
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        downloaded = 0
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):  # 1 MB chunks
                f.write(chunk)
                downloaded += len(chunk)
    print(f" {downloaded / 1e6:.1f} MB")


# ── Tar extraction helper ─────────────────────────────────────────────────────

def _extract_batch(tar_path: str, filenames: set[str], dest_dir: str) -> dict[str, str]:
    """
    Extract all requested files from a tar.xz archive in a single pass.
    Returns a mapping of basename → extracted file path for every file found.
    """
    extracted: dict[str, str] = {}
    with tarfile.open(tar_path, "r:xz") as tar:
        for member in tar:
            name = os.path.basename(member.name)
            if name in filenames:
                tar.extract(member, path=dest_dir)
                extracted[name] = os.path.join(dest_dir, member.name)
    return extracted


# ── CSV initialisation ────────────────────────────────────────────────────────

def _ensure_csv_header() -> None:
    if not os.path.exists(OUT_FILE):
        with open(OUT_FILE, "w") as f:
            f.write(
                "Genome,Species,TAA,TAG,TGA,Total_Valid_Stop_Codons,CDS_Count,"
                "Invalid_Count,TAA_Proportion,TAG_Proportion,TGA_Proportion,"
                "Avg_CDS_Length,GC_Content,Genome_Size_bp,Contig_Count\n"
            )


# ── Main pipeline ─────────────────────────────────────────────────────────────

def run_batch_pipeline() -> None:
    """
    Match FA and bakta samples, discover bakta archive URLs from OSF, then
    download and process one batch at a time, cleaning up after each batch.
    """

    # ── Load list files ───────────────────────────────────────────────────────
    fa_info: dict[str, dict] = {}
    with open(FA_LIST_PATH) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            fa_info[row["sample"]] = row

    bakta_info: dict[str, dict] = {}
    with open(BAKTA_LIST_PATH) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("status") == "PASS":
                bakta_info[row["sample"]] = row

    matched = [s for s in fa_info if s in bakta_info]
    print(f"{len(matched)} samples matched between FA and bakta lists")

    # ── Resume from progress file ─────────────────────────────────────────────
    done = _load_progress()
    todo = [s for s in matched if s not in done]
    print(f"{len(done)} already processed, {len(todo)} remaining\n")
    if not todo:
        print("Nothing to do.")
        return

    # ── Discover bakta archive URLs from OSF ─────────────────────────────────
    print("Fetching bakta archive URLs from OSF ...")
    try:
        bakta_batch_urls = _fetch_bakta_urls_from_osf()
        print(f"  Found: {sorted(bakta_batch_urls)}\n")
    except Exception as e:
        print(f"  OSF lookup failed: {e}")
        print("  Cannot continue without bakta download URLs.\n")
        return

    # ── Group remaining samples by FA batch archive ───────────────────────────
    batches: dict[str, list[str]] = defaultdict(list)
    for sample in todo:
        batches[fa_info[sample]["tar_xz"]].append(sample)

    _ensure_csv_header()
    os.makedirs(TEMP_DIR, exist_ok=True)

    # ── Process one batch at a time ───────────────────────────────────────────
    for fa_tar_name, samples in sorted(batches.items()):
        bakta_tar_name = bakta_info[samples[0]]["tar_xz"]

        print(f"── Batch: {fa_tar_name}  ({len(samples)} samples) ──")

        if bakta_tar_name not in bakta_batch_urls:
            print(f"  [SKIP] {bakta_tar_name} not found in OSF project.\n")
            continue

        fa_tar_path    = os.path.join(TEMP_DIR, fa_tar_name)
        bakta_tar_path = os.path.join(TEMP_DIR, bakta_tar_name)

        # Download both archives for this batch
        try:
            _download(fa_info[samples[0]]["tar_xz_url"], fa_tar_path)
            _download(bakta_batch_urls[bakta_tar_name],  bakta_tar_path)
        except Exception as e:
            print(f"  [ERROR] Download failed: {e}\n")
            for p in [fa_tar_path, bakta_tar_path]:
                if os.path.exists(p):
                    os.remove(p)
            continue

        # Extract all files from both archives in one pass each
        fa_filenames    = {os.path.basename(fa_info[s]["filename_in_tar_xz"]) for s in samples}
        bakta_filenames = {bakta_info[s]["file_name"] for s in samples}

        print(f"  Extracting {len(fa_filenames)} FA files ...")
        fa_extracted    = _extract_batch(fa_tar_path,    fa_filenames,    TEMP_DIR)
        print(f"  Extracting {len(bakta_filenames)} bakta files ...")
        bakta_extracted = _extract_batch(bakta_tar_path, bakta_filenames, TEMP_DIR)

        # Process each sample using the already-extracted files
        for sample in samples:
            fa_filename    = os.path.basename(fa_info[sample]["filename_in_tar_xz"])
            bakta_filename = bakta_info[sample]["file_name"]
            fa_path        = fa_extracted.get(fa_filename)
            bakta_path     = bakta_extracted.get(bakta_filename)

            try:
                if not fa_path or not bakta_path:
                    raise FileNotFoundError(
                        f"Missing extracted file: fa={fa_filename}, bakta={bakta_filename}"
                    )

                c       = count_stop_codons_json(fa_path, bakta_path)
                total   = c["total_valid"]
                avg     = safe_proportion(c["total_cds_length"], c["cds_count"])
                species = fa_info[sample].get("species_sylph", "")

                row = (
                    f"{sample},{species},{c['TAA']},{c['TAG']},{c['TGA']},"
                    f"{total},{c['cds_count']},{c['invalid_count']},"
                    f"{safe_proportion(c['TAA'], total)},"
                    f"{safe_proportion(c['TAG'], total)},"
                    f"{safe_proportion(c['TGA'], total)},"
                    f"{avg},"
                    f"{c['gc_content']},{c['genome_size']},{c['contig_count']}\n"
                )
                with open(OUT_FILE, "a") as f:
                    f.write(row)

                done.add(sample)
                _save_progress(done)
                print(f"  [OK] {sample}  CDS={c['cds_count']}  Valid={total}")

            except Exception as e:
                print(f"  [ERROR] {sample}: {e}")

            finally:
                for p in filter(None, [fa_path, bakta_path]):
                    try:
                        if os.path.exists(p):
                            os.remove(p)
                    except OSError:
                        pass

        # Delete batch archives before starting the next batch
        for p in [fa_tar_path, bakta_tar_path]:
            if os.path.exists(p):
                os.remove(p)
                print(f"  Deleted {os.path.basename(p)}")
        print()

    print(f"Pipeline complete.  Results → {OUT_FILE}")


if __name__ == "__main__":
    run_batch_pipeline()
