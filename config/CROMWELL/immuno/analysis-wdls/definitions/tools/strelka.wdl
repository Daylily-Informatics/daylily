version 1.0

task strelka {
  input {
    File tumor_bam
    File tumor_bam_bai
    File normal_bam
    File normal_bam_bai
    File reference
    File reference_fai
    File reference_dict
    Boolean exome_mode
    Int? cpu_reserved
    # `call_regions` is available to avoid performance issues in specific cases
    # https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md#improving-runtime-for-references-with-many-short-contigs-such-as-grch38
    File? call_regions
    File? call_regions_tbi
  }

  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Float bam_size = size([tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai], "GB")
  Int space_needed_gb = 10 + round(bam_size*2 + reference_size)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    cpu: 4
    docker: "mgibio/strelka-cwl:2.9.9"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    # Ensure bam and bai have same paths, bai timestamp after bam timestamp
    mv ~{tumor_bam} ~{basename(tumor_bam)}; mv ~{tumor_bam_bai} ~{basename(tumor_bam_bai)}
    mv ~{normal_bam} ~{basename(normal_bam)}; mv ~{normal_bam_bai} ~{basename(normal_bam_bai)}

    /usr/bin/perl /usr/bin/docker_helper.pl \
    ~{if defined(cpu_reserved) then cpu_reserved else ""} \
    "$PWD" --tumorBam=~{basename(tumor_bam)} --normalBam=~{basename(normal_bam)} \
    --referenceFasta=~{reference} \
    ~{if defined(call_regions) then "--callRegions=~{call_regions}" else ""} \
    ~{if exome_mode then "--exome" else ""}
  >>>

  output {
    File indels = "results/variants/somatic.indels.vcf.gz"
    File indels_tbi = "results/variants/somatic.indels.vcf.gz.tbi"
    File snvs = "results/variants/somatic.snvs.vcf.gz"
    File snvs_tbi = "results/variants/somatic.snvs.vcf.gz.tbi"
  }
}

workflow wf {
  input {
    File tumor_bam
    File tumor_bam_bai
    File normal_bam
    File normal_bam_bai
    File reference
    File reference_fai
    File reference_dict
    Boolean exome_mode
    Int? cpu_reserved
  }

  call strelka {
    input:
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    exome_mode=exome_mode,
    cpu_reserved=cpu_reserved,
  }
}
