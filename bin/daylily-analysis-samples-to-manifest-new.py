#!/usr/bin/env python3
import os
import sys
import csv
import subprocess
import requests
from pathlib import Path
import boto3
from botocore.exceptions import NoCredentialsError
from collections import defaultdict

def log_info(message):
    print(f"[INFO] {message}")

def log_warn(message):
    print(f"[WARN] {message}")

def log_error(message):
    print(f"[ERROR] {message}")
    exit(1)

def check_file_exists(file_path):
    if file_path.startswith(("http://", "https://")):
        response = requests.head(file_path)
        if response.status_code != 200:
            log_error(f"HTTP file not found: {file_path}")
    elif file_path.startswith("s3://"):
        s3 = boto3.client("s3")
        bucket, key = file_path[5:].split("/", 1)
        try:
            s3.head_object(Bucket=bucket, Key=key)
        except NoCredentialsError:
            log_error("AWS credentials not configured.")
        except Exception as e:
            log_error(f"S3 file not found: {file_path} ({e})")
    else:
        if not os.path.exists(file_path):
            log_error(f"Local file not found: {file_path}")

def determine_sex(n_x, n_y):
    if n_x == 2 and n_y == 0:
        return "female"
    elif n_x == 1 and n_y == 1:
        return "male"
    return "na"

def validate_and_stage_concordance_dir(concordance_dir, stage_target, sample_prefix):
    if concordance_dir == "na" or concordance_dir.startswith("/fsx/data"):
        return concordance_dir
    target_concordance_dir = os.path.join(stage_target, sample_prefix, "concordance_data")
    os.makedirs(target_concordance_dir, exist_ok=True)
    if concordance_dir.startswith(("http://", "https://")):
        subprocess.run(["wget", "-q", "-P", target_concordance_dir, concordance_dir], check=True)
    elif concordance_dir.startswith("s3://"):
        subprocess.run(["aws", "s3", "cp", concordance_dir, target_concordance_dir, "--recursive"], check=True)
    return target_concordance_dir

def validate_subsample_pct(subsample_pct):
    try:
        pct = float(subsample_pct)
        return pct if 0.0 < pct < 1.0 else "na"
    except ValueError:
        return "na"

def generate_analysis_manifest(manifest_file, rows):
    header = [
        "samp", "sample", "sample_lane", "SQ", "RU", "EX", "LANE", "r1_path", "r2_path",
        "biological_sex", "iddna_uid", "concordance_control_path", "is_positive_control",
        "is_negative_control", "sample_type", "merge_single", "external_sample_id",
        "instrument", "lib_prep", "bwa_kmer", "subsample_pct"
    ]
    with open(manifest_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)

def copy_files_to_target(src, dst, link=False):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if src.startswith("s3://"):
        bucket, key = src[5:].split("/", 1)
        boto3.client("s3").download_file(bucket, key, dst)
    elif src.startswith(("http://", "https://")):
        subprocess.run(["wget", "-q", "-O", dst, src], check=True)
    else:
        if link:
            subprocess.run(["ln", "-s", src, dst], check=True)
        else:
            subprocess.run(["cp", src, dst], check=True)

def parse_and_validate_tsv(input_file, stage_target):
    samples = defaultdict(list)
    with open(input_file) as ff:
        next(ff)
        for line in ff:
            cols = line.strip().split("\t")
            key = tuple(cols[:6])
            samples[key].append(cols)

    rows = []
    for sample_key, entries in samples.items():
        is_multi_lane = len(entries) > 1
        lanes = [e[6] for e in entries]
        if is_multi_lane and "0" in lanes:
            log_error(f"Invalid LANE=0 for multi-lane sample: {sample_key}")

        sample_prefix = f"{sample_key[0]}_{sample_key[1]}-{sample_key[2]}_{sample_key[3]}"
        staged_sample_path = os.path.join(stage_target, sample_prefix)
        os.makedirs(staged_sample_path, exist_ok=True)

        if is_multi_lane:
            merged_r1 = os.path.join(staged_sample_path, f"{sample_key[1]}_merged_R1.fastq.gz")
            merged_r2 = os.path.join(staged_sample_path, f"{sample_key[1]}_merged_R2.fastq.gz")
            r1_files, r2_files = zip(*[(e[9], e[10]) for e in entries])

            for f in r1_files + r2_files:
                check_file_exists(f)

            tmp_r1_files = []
            tmp_r2_files = []

            # Download S3 files locally first if they're from S3
            for idx, (r1, r2) in enumerate(zip(r1_files, r2_files)):
                local_r1 = os.path.join(staged_sample_path, f"tmp_{idx}_R1.fastq.gz")
                local_r2 = os.path.join(staged_sample_path, f"tmp_{idx}_R2.fastq.gz")
                copy_files_to_target(r1, local_r1)
                copy_files_to_target(r2, local_r2)
                tmp_r1_files.append(local_r1)
                tmp_r2_files.append(local_r2)

            # Concatenate the downloaded local files
            subprocess.run(f"cat {' '.join(tmp_r1_files)} > {merged_r1}", shell=True, check=True)
            subprocess.run(f"cat {' '.join(tmp_r2_files)} > {merged_r2}", shell=True, check=True)

            # Clean up temporary files
            for tmp_file in tmp_r1_files + tmp_r2_files:
                os.remove(tmp_file)

            rows.append([
                sample_prefix, sample_prefix, sample_prefix, sample_key[3], sample_key[0], sample_key[1],
                "0", merged_r1, merged_r2, determine_sex(int(entries[0][16]), int(entries[0][17])), "na",
                validate_and_stage_concordance_dir(entries[0][8], stage_target, sample_prefix),
                entries[0][14], entries[0][15], sample_key[3], "merge", sample_key[1],
                sample_key[5], sample_key[4], "19", validate_subsample_pct(entries[0][13])
            ])
        else:
            entry = entries[0]
            staged_r1 = os.path.join(staged_sample_path, os.path.basename(entry[9]))
            staged_r2 = os.path.join(staged_sample_path, os.path.basename(entry[10]))
            copy_files_to_target(entry[9], staged_r1, entry[11] == "link_data")
            copy_files_to_target(entry[10], staged_r2, entry[11] == "link_data")

            rows.append([
                sample_prefix, sample_prefix, sample_prefix, sample_key[3], sample_key[0], sample_key[1],
                entry[6], staged_r1, staged_r2, determine_sex(int(entry[16]), int(entry[17])), "na",
                validate_and_stage_concordance_dir(entry[8], stage_target, sample_prefix),
                entry[14], entry[15], sample_key[3], "merge", sample_key[1],
                sample_key[5], sample_key[4], "19", validate_subsample_pct(entry[13])
            ])

    manifest_file = os.path.join(stage_target, "analysis_manifest.csv")
    generate_analysis_manifest(manifest_file, rows)
    log_info(f"Manifest created: {manifest_file}")
    log_info(f"Use this manifest: \n\tcp {manifest_file} config/analysis_manifest.csv")

# Add the following main function to handle command-line arguments and invoke parsing
def main():
    if len(sys.argv) != 3:
        log_error("Usage: script.py <input_tsv> <stage_target>")

    input_file = sys.argv[1]
    stage_target = sys.argv[2]

    parse_and_validate_tsv(input_file, stage_target)

if __name__ == "__main__":
    main()