version 1.0

task germlineCombineVariants {
  input {
    File reference
    File reference_fai
    File reference_dict
    File varscan_vcf
    File varscan_vcf_tbi
    File docm_vcf
    File docm_vcf_tbi
  }

  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Float vcf_size = size([varscan_vcf, varscan_vcf_tbi, docm_vcf, docm_vcf_tbi], "GB")
  Int space_needed_gb = 10 + round(reference_size + 2*vcf_size)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "9GB"
    bootDiskSizeGb: 25
    docker: "mgibio/gatk-cwl:3.6.0"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "combined.vcf.gz"
  command <<<
    /usr/bin/java -Xmx8g -jar /opt/GenomeAnalysisTK.jar -T CombineVariants \
    -genotypeMergeOptions PRIORITIZE \
    --rod_priority_list varscan,docm \
    -o ~{outfile} \
    -R ~{reference} --variant:varscan ~{varscan_vcf} --variant:docm ~{docm_vcf}
  >>>

  output {
    File combined_vcf = outfile
    File combined_vcf_tbi = outfile + ".tbi"
  }
}
