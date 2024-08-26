version 1.0

task pindel {
  input {
    File reference
    File reference_fai
    File reference_dict
    File tumor_bam
    File tumor_bam_bai
    File normal_bam
    File normal_bam_bai
    File region_file
    String tumor_sample_name
    String normal_sample_name
    String? chromosome
    Int insert_size = 400
  }

  Int cores = 4
  Int space_needed_gb = 10 + round(size([reference, reference_fai, reference_dict, normal_bam, normal_bam_bai, tumor_bam, tumor_bam_bai, region_file], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    bootDiskSizeGb: 100
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
    docker: "mgibio/cle:v1.4.2"
    memory: "16GB"
  }

  command <<<
    mv ~{tumor_bam} ~{basename(tumor_bam)}; mv ~{tumor_bam_bai} ~{basename(tumor_bam_bai)}
    mv ~{normal_bam} ~{basename(normal_bam)}; mv ~{normal_bam_bai} ~{basename(normal_bam_bai)}

    echo -e "~{basename(normal_bam)}\t~{insert_size}\t~{normal_sample_name}" > pindel.config
    echo -e "~{basename(tumor_bam)}\t~{insert_size}\t~{tumor_sample_name}" >> pindel.config

    /usr/bin/pindel -i pindel.config -w 30 -T ~{cores} -o all -f ~{reference} \
    ~{if defined(chromosome) then "-c ~{chromosome}" else ""} \
    ~{if defined(region_file) then "-j ~{region_file}" else ""}
  >>>

  output {
    File deletions = "all_D"
    File insertions = "all_SI"
    File tandems = "all_TD"
    File long_insertions = "all_LI"
    File inversions = "all_INV"
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict
    File tumor_bam
    File tumor_bam_bai
    File normal_bam
    File normal_bam_bai
    File region_file
    String tumor_sample_name
    String normal_sample_name
    String? chromosome
    Int insert_size = 400
  }

  call pindel {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    region_file=region_file,
    tumor_sample_name=tumor_sample_name,
    normal_sample_name=normal_sample_name,
    chromosome=chromosome,
    insert_size=insert_size
  }
}
