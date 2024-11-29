#!/usr/bin/env python3

import os
import subprocess
import sys
import shutil
import re

INFO = "\033[1;36m"
WARN = "\033[1;33m"
ERROR = "\033[1;31m"
SUCCESS = "\033[1;32m"
RESET = "\033[0m"

AM_OUT_FILE = "na"


def log_info(message):
    print(f"{INFO}[INFO]{RESET} {message}")


def log_warn(message):
    print(f"{WARN}[WARN]{RESET} {message}")


def log_error(message):
    print(f"{ERROR}[ERROR]{RESET} {message}")
    sys.exit(1)


def log_success(message):
    print(f"{SUCCESS}[SUCCESS]{RESET} {message}")


def check_dependencies():
    required_cmds = ["ssh", "scp"]
    for cmd in required_cmds:
        if not shutil.which(cmd):
            log_error(f"Required command '{cmd}' is not installed.")


def validate_args(mode, input_file):
    if mode not in ["local", "remote"]:
        log_error("First argument must be 'local' or 'remote'.")
    if not os.path.isfile(input_file):
        log_error(f"Sample sheet '{input_file}' does not exist.")


def determine_sex(n_x, n_y):
    if n_x == 2 and n_y == 0:
        return "female"
    elif n_x == 1 and n_y == 1:
        return "male"
    else:
        return "na"


def validate_stage_target(stage_target, mode, cluster_ip=None, pem_file=None, cluster_user=None):
    log_info(f"Validating STAGE_TARGET: {mode} {stage_target}")
    if not stage_target.startswith("/fsx"):
        log_error(f"STAGE_TARGET '{stage_target}' must begin with '/fsx'.")

    if mode == "remote":
        log_info(f"Validating STAGE_TARGET on remote cluster: {stage_target}")
        cmd = f"ssh -i {pem_file} {cluster_user}@{cluster_ip} 'mkdir -p {stage_target}'"
        result = subprocess.run(cmd, shell=True)
        if result.returncode != 0:
            log_error(f"Failed to create STAGE_TARGET directory on remote: {stage_target}")
    else:
        if not os.path.isdir(stage_target):
            log_info(f"Creating STAGE_TARGET directory locally: {stage_target}")
            try:
                os.makedirs(stage_target)
            except Exception as e:
                log_error(f"Failed to create STAGE_TARGET directory: {stage_target} ({e})")


def stage_data_files(r1_src, r2_src, stage_target, sample_prefix, mode, cluster_ip=None, pem_file=None, cluster_user=None):
    r1_dest = os.path.join(stage_target, f"{sample_prefix}.R1.fastq.gz")
    r2_dest = os.path.join(stage_target, f"{sample_prefix}.R2.fastq.gz")

    if mode == "remote":
        log_info(f"Staging R1 to remote: {r1_dest}")
        ssh_cmd = f"ssh -i {pem_file} {cluster_user}@{cluster_ip} 'if [[ -f {r1_dest} ]]; then exit 1; fi'"
        if subprocess.run(ssh_cmd, shell=True).returncode != 0:
            log_error(f"File '{r1_dest}' already exists on remote.")
        scp_cmd = f"scp -i {pem_file} {r1_src} {cluster_user}@{cluster_ip}:{r1_dest}"
        subprocess.run(scp_cmd, shell=True, check=True)
        log_success(f"R1 staged to remote: {r1_dest}")

        log_info(f"Staging R2 to remote: {r2_dest}")
        ssh_cmd = f"ssh -i {pem_file} {cluster_user}@{cluster_ip} 'if [[ -f {r2_dest} ]]; then exit 1; fi'"
        if subprocess.run(ssh_cmd, shell=True).returncode != 0:
            log_error(f"File '{r2_dest}' already exists on remote.")
        scp_cmd = f"scp -i {pem_file} {r2_src} {cluster_user}@{cluster_ip}:{r2_dest}"
        subprocess.run(scp_cmd, shell=True, check=True)
        log_success(f"R2 staged to remote: {r2_dest}")
    else:
        if not os.path.exists(r1_src) or not os.path.exists(r2_src):
            log_error(f"One of the source files does not exist: {r1_src}, {r2_src}")
        shutil.copy(r1_src, r1_dest)
        log_success(f"R1 staged locally: {r1_dest}")
        shutil.copy(r2_src, r2_dest)
        log_success(f"R2 staged locally: {r2_dest}")


def parse_and_validate_tsv(input_file, mode, cluster_ip=None, pem_file=None, cluster_user=None):
    global AM_OUT_FILE
    with open(input_file, "r") as f:
        lines = f.readlines()

    for line in lines[1:]:
        cols = line.strip().split("\t")
        if len(cols) < 17:
            log_error(f"Invalid line in TSV: {line}")
        (RUN_ID, SAMPLE_ID, SAMPLE_TYPE, LIB_PREP, SEQ_PLATFORM, LANE, SEQBC_ID, 
         PATH_TO_CONCORDANCE_DATA_DIR, R1_FQ, R2_FQ, STAGE_DIRECTIVE, STAGE_TARGET, 
         SUBSAMPLE_PCT, IS_POS_CTRL, IS_NEG_CTRL, N_X, N_Y) = cols[:17]

        if LANE != "0":
            log_error(f"Invalid LANE value '{LANE}'. Only '0' is allowed.")
        if IS_POS_CTRL != "na" or IS_NEG_CTRL != "na":
            log_error("Positive or Negative Control columns are not enabled yet.")

        log_info(f"Processing: {line.strip()}")
        validate_stage_target(STAGE_TARGET, mode, cluster_ip, pem_file, cluster_user)

        if STAGE_DIRECTIVE == "stage_data":
            stage_data_files(R1_FQ, R2_FQ, STAGE_TARGET, f"{RUN_ID}_{SAMPLE_ID}_{SAMPLE_TYPE}_{LANE}-{SEQBC_ID}", mode, cluster_ip, pem_file, cluster_user)


def main():
    if len(sys.argv) < 3:
        log_error("Usage: <local|remote> <input_file> [<cluster_ip> <pem_file> <cluster_user>]")

    mode = sys.argv[1]
    input_file = sys.argv[2]
    cluster_ip = sys.argv[3] if len(sys.argv) > 3 else None
    pem_file = sys.argv[4] if len(sys.argv) > 4 else None
    cluster_user = sys.argv[5] if len(sys.argv) > 5 else None

    check_dependencies()
    validate_args(mode, input_file)
    parse_and_validate_tsv(input_file, mode, cluster_ip, pem_file, cluster_user)
    log_success(f"Manifest generated: {AM_OUT_FILE}")


if __name__ == "__main__":
    main()
