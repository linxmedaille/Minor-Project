import csv
import json
import os
import tarfile

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
    tar_files_by_name = {}
    pending_urls = [storage_url]

    while len(pending_urls) > 0:
        url = pending_urls.pop(0)
        response = requests.get(url, params={"page[size]": 100}, timeout=30)
        response.raise_for_status()
        payload = response.json()

        for entry in payload["data"]:
            name = entry["attributes"]["name"]
            entry_type = entry["attributes"]["kind"]

            if entry_type == "file" and name.endswith(".tar.xz"):
                tar_files_by_name[name] = entry["links"]["download"]
            elif entry_type == "folder":
                pending_urls.append(entry["relationships"]["files"]["links"]["related"]["href"])

        next_page_url = payload["links"].get("next")
        if next_page_url:
            pending_urls.append(next_page_url)

    return tar_files_by_name


def get_bakta_urls():
    storage_url = f"{OSF_API}/nodes/{BAKTA_OSF_NODE}/files/osfstorage/"
    all_files = list_osf_files(storage_url)

    bakta_files = {}
    for name in all_files:
        if "bakta" in name:
            bakta_files[name] = all_files[name]
    return bakta_files


def load_progress():
    if not os.path.exists(PROGRESS_FILE):
        return []

    with open(PROGRESS_FILE) as progress_file:
        content = progress_file.read().strip()

    if content == "":
        return []

    return json.loads(content)


def save_progress(done_list):
    with open(PROGRESS_FILE, "w") as progress_file:
        json.dump(sorted(done_list), progress_file, indent=2)


def download_file(url, dest):
    name = os.path.basename(dest)
    print(f"  Downloading {name} ...", end="", flush=True)

    one_megabyte = 1024 * 1024
    downloaded = 0

    with requests.get(url, stream=True, timeout=60) as response:
        response.raise_for_status()
        with open(dest, "wb") as output_file:
            for chunk in response.iter_content(chunk_size=one_megabyte):
                output_file.write(chunk)
                downloaded += len(chunk)

    print(f" {round(downloaded / 1000000, 1)} MB")


def extract_files(tar_path, filenames_wanted, dest_dir):
    extracted = {}

    with tarfile.open(tar_path, "r:xz") as tar:
        for member in tar:
            name = os.path.basename(member.name)
            if name in filenames_wanted:
                tar.extract(member, path=dest_dir)
                extracted[name] = os.path.join(dest_dir, member.name)

    return extracted


def init_csv():
    if os.path.exists(OUT_FILE):
        return

    with open(OUT_FILE, "w") as csv_file:
        csv_file.write(
            "Genome,Species,TAA,TAG,TGA,Total_Valid_Stop_Codons,CDS_Count,"
            "Invalid_Count,TAA_Proportion,TAG_Proportion,TGA_Proportion,"
            "Avg_CDS_Length,GC_Content,Genome_Size_bp,Contig_Count\n"
        )


def group_samples_by_tar(sample_list, fa_info):
    batches = {}
    for sample in sample_list:
        tar_name = fa_info[sample]["tar_xz"]
        if tar_name not in batches:
            batches[tar_name] = []
        batches[tar_name].append(sample)
    return batches


def unique_names(name_list):
    result = []
    for name in name_list:
        if name not in result:
            result.append(name)
    return result


def delete_if_exists(path):
    if not path:
        return
    try:
        if os.path.exists(path):
            os.remove(path)
    except OSError:
        pass


def run_pipeline():
    fa_info = {}
    with open(FA_LIST_PATH) as fa_list_file:
        for row in csv.DictReader(fa_list_file, delimiter="\t"):
            fa_info[row["sample"]] = row

    bakta_info = {}
    with open(BAKTA_LIST_PATH) as bakta_list_file:
        for row in csv.DictReader(bakta_list_file, delimiter="\t"):
            if row.get("status") == "PASS":
                bakta_info[row["sample"]] = row

    matched = []
    for sample in fa_info:
        if sample in bakta_info:
            matched.append(sample)
    print(f"{len(matched)} samples matched between FA and bakta lists")

    done = load_progress()
    todo = []
    for sample in matched:
        if sample not in done:
            todo.append(sample)

    print(f"{len(done)} already processed, {len(todo)} remaining")
    if len(todo) == 0:
        print("Nothing to do.")
        return

    print("Fetching bakta archive URLs from OSF ...")
    try:
        bakta_batch_urls = get_bakta_urls()
        print(f"  Found: {sorted(bakta_batch_urls)}")
    except Exception as error:
        print(f"  OSF lookup failed: {error}")
        print("  Cannot continue without bakta download URLs.")
        return

    batches = group_samples_by_tar(todo, fa_info)

    init_csv()
    os.makedirs(TEMP_DIR, exist_ok=True)

    for fa_tar_name in sorted(batches):
        samples = batches[fa_tar_name]
        bakta_tar_name = bakta_info[samples[0]]["tar_xz"]

        print(f"Batch: {fa_tar_name}  ({len(samples)} samples)")

        if bakta_tar_name not in bakta_batch_urls:
            print(f"  [SKIP] {bakta_tar_name} not found in OSF project.")
            continue

        fa_tar_path = os.path.join(TEMP_DIR, fa_tar_name)
        bakta_tar_path = os.path.join(TEMP_DIR, bakta_tar_name)

        try:
            download_file(fa_info[samples[0]]["tar_xz_url"], fa_tar_path)
            download_file(bakta_batch_urls[bakta_tar_name], bakta_tar_path)
        except Exception as error:
            print(f"  [ERROR] Download failed: {error}")
            delete_if_exists(fa_tar_path)
            delete_if_exists(bakta_tar_path)
            continue

        fa_filenames = []
        bakta_filenames = []
        for sample in samples:
            fa_filenames.append(os.path.basename(fa_info[sample]["filename_in_tar_xz"]))
            bakta_filenames.append(bakta_info[sample]["file_name"])
        fa_filenames = unique_names(fa_filenames)
        bakta_filenames = unique_names(bakta_filenames)

        print(f"  Extracting {len(fa_filenames)} FA files ...")
        fa_extracted = extract_files(fa_tar_path, fa_filenames, TEMP_DIR)
        print(f"  Extracting {len(bakta_filenames)} bakta files ...")
        bakta_extracted = extract_files(bakta_tar_path, bakta_filenames, TEMP_DIR)

        for sample in samples:
            fa_filename = os.path.basename(fa_info[sample]["filename_in_tar_xz"])
            bakta_filename = bakta_info[sample]["file_name"]
            fa_path = fa_extracted.get(fa_filename)
            bakta_path = bakta_extracted.get(bakta_filename)

            if not fa_path or not bakta_path:
                print(f"  [ERROR] {sample}: missing extracted file "
                      f"(fa={fa_filename}, bakta={bakta_filename})")
            else:
                try:
                    counts = count_stop_codons_json(fa_path, bakta_path)
                    total = counts["total_valid"]
                    avg_length = safe_proportion(counts["total_cds_length"], counts["cds_count"])
                    species = fa_info[sample].get("species_sylph", "")

                    values = [
                        sample,
                        species,
                        str(counts["TAA"]),
                        str(counts["TAG"]),
                        str(counts["TGA"]),
                        str(total),
                        str(counts["cds_count"]),
                        str(counts["invalid_count"]),
                        safe_proportion(counts["TAA"], total),
                        safe_proportion(counts["TAG"], total),
                        safe_proportion(counts["TGA"], total),
                        avg_length,
                        str(counts["gc_content"]),
                        str(counts["genome_size"]),
                        str(counts["contig_count"]),
                    ]

                    with open(OUT_FILE, "a") as csv_file:
                        csv_file.write(",".join(values) + "\n")

                    done.append(sample)
                    save_progress(done)
                    print(f"  [OK] {sample}  CDS={counts['cds_count']}  Valid={total}")

                except Exception as error:
                    print(f"  [ERROR] {sample}: {error}")

            delete_if_exists(fa_path)
            delete_if_exists(bakta_path)

        for tar_path in [fa_tar_path, bakta_tar_path]:
            if os.path.exists(tar_path):
                os.remove(tar_path)
                print(f"  Deleted {os.path.basename(tar_path)}")

    print(f"Pipeline complete.  Results -> {OUT_FILE}")


if __name__ == "__main__":
    run_pipeline()
