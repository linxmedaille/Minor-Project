import csv
import json
import os
import tarfile
from collections import defaultdict

import requests

from analysis import count_stop_codons_json, safe_proportion, OUT_FILE


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FA_LIST_PATH = os.path.join(BASE_DIR, "DataFiles", "Lists", "file_list.r0.2.v2.tsv")
BAKTA_LIST_PATH = os.path.join(BASE_DIR, "DataFiles", "Lists", "atb.bakta.r0.2.status.tsv")
PROGRESS_FILE = os.path.join(BASE_DIR, "pipeline_progress.json")
TEMP_DIR = os.path.join(BASE_DIR, "temp_downloads")

OSF_API = "https://api.osf.io/v2"
BAKTA_OSF_NODE = "fqwgr"


def list_osf_files(storage_url):
    found = {}
    url = storage_url

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
                found.update(list_osf_files(sub_url))

        url = payload["links"].get("next")

    return found


def get_bakta_urls():
    storage_url = f"{OSF_API}/nodes/{BAKTA_OSF_NODE}/files/osfstorage/"
    all_files   = list_osf_files(storage_url)
    return {name: url for name, url in all_files.items() if "bakta" in name}


def load_progress():
    if os.path.exists(PROGRESS_FILE):
        with open(PROGRESS_FILE) as f:
            content = f.read().strip()
            if content:
                return set(json.loads(content))
    return set()


def save_progress(done):
    with open(PROGRESS_FILE, "w") as f:
        json.dump(sorted(done), f, indent=2)


def download_file(url, dest):
    name = os.path.basename(dest)
    print(f"  Downloading {name} ...", end="", flush=True)
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        downloaded = 0
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
                downloaded += len(chunk)
    print(f" {downloaded / 1e6:.1f} MB")


def extract_files(tar_path, filenames, dest_dir):
    extracted = {}
    with tarfile.open(tar_path, "r:xz") as tar:
        for member in tar:
            name = os.path.basename(member.name)
            if name in filenames:
                tar.extract(member, path=dest_dir)
                extracted[name] = os.path.join(dest_dir, member.name)
    return extracted


def init_csv():
    if not os.path.exists(OUT_FILE):
        with open(OUT_FILE, "w") as f:
            f.write(
                "Genome,Species,TAA,TAG,TGA,Total_Valid_Stop_Codons,CDS_Count,"
                "Invalid_Count,TAA_Proportion,TAG_Proportion,TGA_Proportion,"
                "Avg_CDS_Length,GC_Content,Genome_Size_bp,Contig_Count\n"
            )


def run_pipeline():
    fa_info = {}
    with open(FA_LIST_PATH) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            fa_info[row["sample"]] = row

    bakta_info = {}
    with open(BAKTA_LIST_PATH) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("status") == "PASS":
                bakta_info[row["sample"]] = row

    matched = [s for s in fa_info if s in bakta_info]
    print(f"{len(matched)} samples matched between FA and bakta lists")

    done = load_progress()
    todo = [s for s in matched if s not in done]
    print(f"{len(done)} already processed, {len(todo)} remaining\n")
    if not todo:
        print("Nothing to do.")
        return

    print("Fetching bakta archive URLs from OSF ...")
    try:
        bakta_batch_urls = get_bakta_urls()
        print(f"  Found: {sorted(bakta_batch_urls)}\n")
    except Exception as e:
        print(f"  OSF lookup failed: {e}\n  Cannot continue without bakta download URLs.\n")
        return

    batches = defaultdict(list)
    for sample in todo:
        batches[fa_info[sample]["tar_xz"]].append(sample)

    init_csv()
    os.makedirs(TEMP_DIR, exist_ok=True)

    for fa_tar_name, samples in sorted(batches.items()):
        bakta_tar_name = bakta_info[samples[0]]["tar_xz"]

        print(f"Batch: {fa_tar_name}  ({len(samples)} samples)")

        if bakta_tar_name not in bakta_batch_urls:
            print(f"  [SKIP] {bakta_tar_name} not found in OSF project.\n")
            continue

        fa_tar_path = os.path.join(TEMP_DIR, fa_tar_name)
        bakta_tar_path = os.path.join(TEMP_DIR, bakta_tar_name)

        try:
            download_file(fa_info[samples[0]]["tar_xz_url"], fa_tar_path)
            download_file(bakta_batch_urls[bakta_tar_name],  bakta_tar_path)
        except Exception as e:
            print(f"  [ERROR] Download failed: {e}\n")
            for p in [fa_tar_path, bakta_tar_path]:
                if os.path.exists(p):
                    os.remove(p)
            continue

        fa_filenames = {os.path.basename(fa_info[s]["filename_in_tar_xz"]) for s in samples}
        bakta_filenames = {bakta_info[s]["file_name"] for s in samples}

        print(f"  Extracting {len(fa_filenames)} FA files ...")
        fa_extracted = extract_files(fa_tar_path, fa_filenames, TEMP_DIR)
        print(f"  Extracting {len(bakta_filenames)} bakta files ...")
        bakta_extracted = extract_files(bakta_tar_path, bakta_filenames, TEMP_DIR)

        for sample in samples:
            fa_filename = os.path.basename(fa_info[sample]["filename_in_tar_xz"])
            bakta_filename = bakta_info[sample]["file_name"]
            fa_path = fa_extracted.get(fa_filename)
            bakta_path = bakta_extracted.get(bakta_filename)

            try:
                if not fa_path or not bakta_path:
                    raise FileNotFoundError(
                        f"Missing extracted file: fa={fa_filename}, bakta={bakta_filename}"
                    )

                c = count_stop_codons_json(fa_path, bakta_path)
                total = c["total_valid"]
                avg = safe_proportion(c["total_cds_length"], c["cds_count"])
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
                save_progress(done)
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

        for p in [fa_tar_path, bakta_tar_path]:
            if os.path.exists(p):
                os.remove(p)
                print(f"  Deleted {os.path.basename(p)}")
        print()

    print(f"Pipeline complete.  Results → {OUT_FILE}")


if __name__ == "__main__":
    run_pipeline()
