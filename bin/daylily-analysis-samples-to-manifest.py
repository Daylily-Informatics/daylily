#!/usr/bin/env python3
import os
import csv
import subprocess
import requests
from pathlib import Path
import boto3
from botocore.exceptions import NoCredentialsError

def log_info(message):
    print(f"[INFO] {message}")

def log_warn(message):
    print(f"[WARN] {message}")

def log_error(message):
    print(f"[ERROR] {message}")
    exit(1)

def check_file_exists(file_path):
    if file_path.startswith("http://") or file_path.startswith("https://"):
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


def validate_and_stage_concordance_dir(concordance_dir, stage_target, sample_prefix, aws_profile=None):
    if concordance_dir == "na" or concordance_dir.startswith("/fsx/data"):
        return concordance_dir
    stage_path = os.path.join(stage_target, f"{sample_prefix}")
    os.makedirs(stage_path, exist_ok=True)
    target_concordance_dir = os.path.join(stage_path, "concordance_data")
    os.makedirs(target_concordance_dir, exist_ok=True)
    if concordance_dir.startswith("http://") or concordance_dir.startswith("https://"):
        log_info(f"Downloading concordance data from {concordance_dir}")
        subprocess.run(["wget", "-q", "-P", target_concordance_dir, concordance_dir], check=True)
    elif concordance_dir.startswith("s3://"):
        log_info(f"Downloading concordance data from S3: {concordance_dir}")
        subprocess.run(["aws", "s3", "cp", concordance_dir, target_concordance_dir, "--profile", aws_profile, "--recursive"], check=True)
    return target_concordance_dir

def validate_subsample_pct(subsample_pct):
    try:
        subsample_pct = float(subsample_pct)
        if subsample_pct <= 0.0 or subsample_pct > 1.0:
            return "na"
        return subsample_pct if subsample_pct < 1.0 else "na"
    except ValueError:
        return "na"

def create_remote_directory(cluster_ip, pem_file, cluster_user, remote_path):
    command = f"ssh -i {pem_file} {cluster_user}@{cluster_ip} 'mkdir -p {remote_path}'"
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        log_error(f"Failed to create remote directory: {remote_path}")

def validate_path(mode, cluster_ip, pem_file, cluster_user, path):
    if mode == "remote":
        create_remote_directory(cluster_ip, pem_file, cluster_user, path)
    else:
        os.makedirs(path, exist_ok=True)

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


def copy_files_to_target(file_path, target_path, mode, cluster_ip=None, pem_file=None, cluster_user=None, link=None):
    """Copy a file to the target path, handling local and remote scenarios."""

    try:
        if mode == "remote":
            remote_command = f"scp -i {pem_file} {file_path} {cluster_user}@{cluster_ip}:{target_path}"
            subprocess.run(remote_command, shell=True, check=True)
        else:
            os.makedirs(os.path.dirname(target_path), exist_ok=True)
            if file_path.startswith("http://") or file_path.startswith("https://"):
                subprocess.run(["wget", "-q", "-O", target_path, file_path], check=True)
                log_info(f"Downloaded file {file_path} to {target_path}")
            elif file_path.startswith("s3://"):
                bucket, key = file_path[5:].split("/", 1)
                s3 = boto3.client("s3")
                s3.download_file(bucket, key, target_path)
                log_info(f"Downloaded file {file_path} to {target_path}")
            else:
                if link == 'link':
                    subprocess.run(["ln", "-s", file_path, target_path], check=True)
                else:
                    subprocess.run(["cp", file_path, target_path], check=True)
                    
        log_info(f"Copied file {file_path} to {target_path}")
    except Exception as e:
        log_error(f"Error copying file {file_path} to {target_path}: {e}")

def parse_and_validate_tsv(input_file, mode, cluster_ip=None, pem_file=None, cluster_user=None, aws_profile=None):
    with open(input_file, "r") as ff:
        linesf = ff.readlines()

    stage_target = "na"
    runn = "na"
    rows = []
    stage_directive = ""
    for line in linesf[1:]:
        cols = line.strip().split("\t")
        if len(cols) < 17:
            log_error(f"Invalid line in TSV: {line}")
        (RUN_ID, SAMPLE_ID, SAMPLE_ANNO, SAMPLE_TYPE, LIB_PREP, SEQ_PLATFORM, LANE, SEQBC_ID,
         PATH_TO_CONCORDANCE_DATA_DIR, R1_FQ, R2_FQ, STAGE_DIRECTIVE, STAGE_TARGET,
         SUBSAMPLE_PCT, IS_POS_CTRL, IS_NEG_CTRL, N_X, N_Y) = cols[:18]
        runn = RUN_ID
        stage_target = STAGE_TARGET
        stage_directive = STAGE_DIRECTIVE

    # Ensure stage_target exists
    if stage_target == "is_staged":
        if mode == "remote":
            raise Exception("mode==REMOTE is not supported with 'is_staged' data") #could be however
        else:
            if not os.path.exists(stage_target):
                log_error(f"Stage target directory does not exist: {stage_target}")
                raise Exception("Stage target directory does not exist")
    else:
        try:
            if mode == "remote":
                create_remote_directory(cluster_ip, pem_file, cluster_user, stage_target)
            else:
                os.makedirs(stage_target, exist_ok=True)
        except Exception as e:
            error_file = os.path.join(stage_target, f"{runn}_error.txt")
            with open(error_file, "w") as ef:
                ef.write(f"Error creating stage target: {e}\n")
            log_error(f"Error creating stage target directory: {e}")

    try:
        for line in linesf[1:]:
            cols = line.strip().split("\t")
            (RUN_ID, SAMPLE_ID, SAMPLE_ANNO, SAMPLE_TYPE, LIB_PREP, SEQ_PLATFORM, LANE, SEQBC_ID,
            PATH_TO_CONCORDANCE_DATA_DIR, R1_FQ, R2_FQ, STAGE_DIRECTIVE, STAGE_TARGET,
            SUBSAMPLE_PCT, IS_POS_CTRL, IS_NEG_CTRL, N_X, N_Y) = cols[:18]

            sample_prefix = f"{RUN_ID}_{SAMPLE_ID}-{SAMPLE_ANNO}_{SAMPLE_TYPE}_{LANE}-{SEQBC_ID}"
            staged_sample_path = os.path.join(stage_target, sample_prefix)


            # Validate and stage files
            check_file_exists(R1_FQ)
            check_file_exists(R2_FQ)
            
            # Ensure the directory exists
            if mode == "remote":
                create_remote_directory(cluster_ip, pem_file, cluster_user, staged_sample_path)
            else:
                os.makedirs(staged_sample_path, exist_ok=True)

            # Copy FASTQ files
            staged_r1 = os.path.join(staged_sample_path, os.path.basename(R1_FQ))
            staged_r2 = os.path.join(staged_sample_path, os.path.basename(R2_FQ))
            if STAGE_DIRECTIVE == "stage_data":
                copy_files_to_target(R1_FQ, staged_r1, mode, cluster_ip, pem_file, cluster_user)
                copy_files_to_target(R2_FQ, staged_r2, mode, cluster_ip, pem_file, cluster_user)
            
                # Copy analysis_samples.tsv to stage_target
                analysis_samples_target = os.path.join(stage_target, f"{runn}_analysis_samples.tsv")
                copy_files_to_target(input_file, analysis_samples_target, mode, cluster_ip, pem_file, cluster_user)
            elif STAGE_DIRECTIVE == "link_data":
                copy_files_to_target(R1_FQ, staged_r1, mode, cluster_ip, pem_file, cluster_user, 'link')
                copy_files_to_target(R2_FQ, staged_r2, mode, cluster_ip, pem_file, cluster_user, 'link')

                # Copy analysis_samples.tsv to stage_target
                analysis_samples_target = os.path.join(stage_target, f"{runn}_analysis_samples.tsv")
                copy_files_to_target(input_file, analysis_samples_target, mode, cluster_ip, pem_file, cluster_user)

            else:
                fileex=True
                if not os.path.exists(R1_FQ):
                    log_error(f"FASTQ file not found: {R1_FQ}")
                    fileex=False
                if not os.path.exists(R2_FQ):
                    log_error(f"FASTQ file not found: {R2_FQ}")
                    fileex=False    
                if not fileex:
                    log_error(f"FASTQ file not found")
                    raise Exception("One of fastqs not found", R1_FQ, R2_FQ)


            concordance_dir = validate_and_stage_concordance_dir(PATH_TO_CONCORDANCE_DATA_DIR, staged_sample_path, sample_prefix)
            subsample_pct = validate_subsample_pct(SUBSAMPLE_PCT)
            biological_sex = determine_sex(int(N_X), int(N_Y))

            row = [
                sample_prefix, sample_prefix, sample_prefix, SAMPLE_TYPE, RUN_ID, SAMPLE_ID,
                f"{LANE}-{SEQBC_ID}", staged_r1, staged_r2, biological_sex, "na",
                concordance_dir, IS_POS_CTRL, IS_NEG_CTRL, SAMPLE_TYPE, "merge",
                SAMPLE_ID, SEQ_PLATFORM, LIB_PREP, "19", subsample_pct
            ]
            rows.append(row)

        # Generate analysis manifest
        manifest_file_tmp = os.path.join('./', f"{runn}_analysis_manifest.csv")
        manifest_file = os.path.join(stage_target, f"{runn}_analysis_manifest.csv")

        if os.path.exists(manifest_file) or os.path.exists(manifest_file_tmp):
            log_error(f"Manifest file already exists: {manifest_file}")
            raise Exception("Manifest file already exists")
        generate_analysis_manifest(manifest_file_tmp, rows)

        if STAGE_DIRECTIVE == "stage_data":
            copy_files_to_target(manifest_file_tmp, manifest_file, mode, cluster_ip, pem_file, cluster_user)
            

        # Create success sentinel
        success_file = os.path.join(stage_target, f"{runn}_success.txt")
        if mode == "remote":
            remote_command = f"ssh -i {pem_file} {cluster_user}@{cluster_ip} 'echo \"Process completed successfully\" > {success_file}'"
            subprocess.run(remote_command, shell=True, check=True)
        else:
            with open(success_file, "w") as sf:
                sf.write("Process completed successfully\n")
        log_info(f"Success file created: {success_file}")
        log_info(f"Manifest generated: {manifest_file_tmp} & {manifest_file}")
        log_info(f"\nTo use the manifest just created, copy it: cp  {manifest_file_tmp} config/analysis_manifest.csv\n")

    except Exception as e:
        # Create error sentinel
        error_file = os.path.join(stage_target, f"{runn}_error.txt")
        if mode == "remote":
            remote_command = f"ssh -i {pem_file} {cluster_user}@{cluster_ip} 'echo \"Error during processing: {e}\" > {error_file}'"
            subprocess.run(remote_command, shell=True, check=True)
        else:
            with open(error_file, "w") as ef:
                ef.write(f"Error during processing: {e}\n")
        log_error(f"Error processing the TSV file: {e}")

def main():
    import sys
    if len(sys.argv) < 3:
        log_error("Usage: python3 script.py local|remote input_file aws_profile [cluster_ip] [pem_file] [cluster_user]")
    mode = sys.argv[1]
    input_file = sys.argv[2]
    aws_profile = sys.argv[3] 
    os.environ["AWS_PROFILE"] = aws_profile  # Set AWS profile
    cluster_ip = sys.argv[4] if mode == "remote" else None
    pem_file = sys.argv[5] if mode == "remote" else None
    cluster_user = sys.argv[6] if mode == "remote" else None

    parse_and_validate_tsv(input_file, mode, cluster_ip, pem_file, cluster_user, aws_profile)


if __name__ == "__main__":
    main()


