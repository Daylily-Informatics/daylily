version 1.0

task bamReadcount {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    String sample
    File vcf
    Int min_mapping_quality = 0
    Int min_base_quality = 20
    String prefix = "NOPREFIX"
  }

  Int space_needed_gb = 10 + round(size([bam, bam_bai, reference, reference_fai, reference_dict, vcf], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/bam_readcount_helper-cwl:1.1.1"
    memory: "16GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String stdout_file = sample + "_bam_readcount.tsv"
  String prefixed_sample = (if prefix == "NOPREFIX" then "" else (prefix + "_")) + sample
  command <<<
    mv ~{bam} ~{basename(bam)}; mv ~{bam_bai} ~{basename(bam_bai)}

    /usr/bin/python -c '
    import sys
    import os
    from cyvcf2 import VCF
    import tempfile
    import csv
    from subprocess import Popen, PIPE

    def generate_region_list(hash):
        fh = tempfile.NamedTemporaryFile("w", delete=False)
        writer = csv.writer(fh, delimiter="\t")
        for chr, positions in hash.items():
            for pos in sorted(positions.keys()):
                writer.writerow([chr, pos, pos])
        fh.close()
        return fh.name

    def filter_sites_in_hash(region_list, bam_file, ref_fasta, prefixed_sample, output_dir, insertion_centric, map_qual, base_qual):
        bam_readcount_cmd = ["/usr/bin/bam-readcount", "-f", ref_fasta, "-l", region_list, "-w", "0", "-b", str(base_qual), "-q", str(map_qual)]
        if insertion_centric:
            bam_readcount_cmd.append("-i")
            output_file = os.path.join(output_dir, prefixed_sample + "_bam_readcount_indel.tsv")
        else:
            output_file = os.path.join(output_dir, prefixed_sample + "_bam_readcount_snv.tsv")
        bam_readcount_cmd.append(bam_file)
        execution = Popen(bam_readcount_cmd, stdout=PIPE, stderr=PIPE)
        stdout, stderr = execution.communicate()
        if execution.returncode == 0:
            with open(output_file, "wb") as output_fh:
                output_fh.write(stdout)
        else:
            sys.exit(stderr)

    min_base_qual = ~{min_base_quality}
    min_mapping_qual = ~{min_mapping_quality}
    output_dir = os.environ["PWD"]
    sample = "~{sample}"
    ref_fasta = "~{reference}"
    bam_file = "~{basename(bam)}"
    prefixed_sample = "~{prefixed_sample}"
    vcf_filename = "~{vcf}"

    vcf_file = VCF(vcf_filename)
    sample_index = vcf_file.samples.index(sample)

    rc_for_indel = {}
    rc_for_snp   = {}
    for variant in vcf_file:
        ref = variant.REF
        chr = variant.CHROM
        start = variant.start
        end = variant.end
        pos = variant.POS
        for var in variant.ALT:
            if len(ref) > 1 or len(var) > 1:
                #it is an indel or mnp
                if len(ref) == len(var) or (len(ref) > 1 and len(var) > 1):
                    sys.stderr.write("Complex variant or MNP will be skipped: %s\t%s\t%s\t%s\n" % (chr, pos, ref , var))
                    continue
                elif len(ref) > len(var):
                    #it is a deletion
                    pos += 1
                    unmodified_ref = ref
                    ref = unmodified_ref[1]
                    var = "-%s" % unmodified_ref[1:]
                else:
                    #it is an insertion
                    var = "+%s" % var[1:]
                if chr not in rc_for_indel:
                    rc_for_indel[chr] = {}
                if pos not in rc_for_indel[chr]:
                    rc_for_indel[chr][pos] = {}
                if ref not in rc_for_indel[chr][pos]:
                    rc_for_indel[chr][pos][ref] = {}
                rc_for_indel[chr][pos][ref] = variant
            else:
                #it is a SNP
                if chr not in rc_for_snp:
                    rc_for_snp[chr] = {}
                if pos not in rc_for_snp[chr]:
                    rc_for_snp[chr][pos] = {}
                if ref not in rc_for_snp[chr][pos]:
                    rc_for_snp[chr][pos][ref] = {}
                rc_for_snp[chr][pos][ref] = variant

    if len(rc_for_snp.keys()) > 0:
        region_file = generate_region_list(rc_for_snp)
        filter_sites_in_hash(region_file, bam_file, ref_fasta, prefixed_sample, output_dir, False, min_mapping_qual, min_base_qual)
    else:
        output_file = os.path.join(output_dir, prefixed_sample + "_bam_readcount_snv.tsv")
        open(output_file, "w").close()

    if len(rc_for_indel.keys()) > 0:
        region_file = generate_region_list(rc_for_indel)
        filter_sites_in_hash(region_file, bam_file, ref_fasta, prefixed_sample, output_dir, True, min_mapping_qual, min_base_qual)
    else:
        output_file = os.path.join(output_dir, prefixed_sample + "_bam_readcount_indel.tsv")
        open(output_file, "w").close()
    ' > ~{stdout_file}
  >>>

  output {
    File snv_bam_readcount_tsv = prefixed_sample + "_bam_readcount_snv.tsv"
    File indel_bam_readcount_tsv = prefixed_sample + "_bam_readcount_indel.tsv"
  }
}

workflow wf {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    String sample
    File vcf
    Int? min_mapping_quality
    Int? min_base_quality
    String? prefix
  }

  call bamReadcount {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    sample=sample,
    vcf=vcf,
    min_mapping_quality=min_mapping_quality,
    min_base_quality=min_base_quality,
    prefix=prefix
  }
}
