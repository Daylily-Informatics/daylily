version 1.0

task readBackedPhasing {
  input {
    File bam
    File bam_index
    File reference
    File reference_fai
    File reference_dict
    File vcf
    File vcf_tbi
  }

  Int space_needed_gb = 10 + round(size([bam, bam_index, reference, reference_fai, reference_dict, vcf, vcf_tbi], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/gatk-cwl:3.6.0"
    memory: "9GB"
    bootDiskSizeGb: 25
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "phased.vcf"
  command <<<
    echo "~{bam_index}"
    head -n 1 ~{bam_index}
    /usr/bin/java -Xmx8g -jar /opt/GenomeAnalysisTK.jar -T ReadBackedPhasing \
    -L ~{vcf} -o ~{outfile} \
    -R ~{reference} \
    -I ~{bam} \
    -V ~{vcf}
  >>>

  output {
    File phased_vcf = outfile
  }
}

workflow wf {
  input {
    File bam
    File bam_index
    File reference
    File reference_fai
    File reference_dict
    File vcf
    File vcf_tbi
  }
  call readBackedPhasing {
    input:
    bam=bam,
    bam_index=bam_index,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    vcf=vcf,
    vcf_tbi=vcf_tbi
  }
}
