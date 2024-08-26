version 1.0

task combineVariants {
  input {
    File reference
    File reference_fai
    File reference_dict
    File mutect_vcf
    File mutect_vcf_tbi
    File varscan_vcf
    File varscan_vcf_tbi
    File strelka_vcf
    File strelka_vcf_tbi
  }

  Float ref_size = size([reference, reference_fai, reference_dict], "GB")
  Float mutect_size = size([mutect_vcf, mutect_vcf_tbi], "GB")
  Float varscan_size = size([varscan_vcf, varscan_vcf_tbi], "GB")
  Float strelka_size = size([strelka_vcf, strelka_vcf_tbi], "GB")
  Int space_needed_gb = 10 + round(ref_size + mutect_size + varscan_size + strelka_size)*2
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "9GB"
    docker:  "mgibio/gatk-cwl:3.6.0"
    bootDiskSizeGb: 25
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "combined.vcf.gz"
  command <<<
    /usr/bin/java -Xmx8g -jar /opt/GenomeAnalysisTK.jar -T CombineVariants \
    -genotypeMergeOptions PRIORITIZE \
    --rod_priority_list mutect,varscan,strelka \
    -o ~{outfile} \
    -R ~{reference} \
    --variant:mutect ~{mutect_vcf} \
    --variant:varscan ~{varscan_vcf} \
    --variant:strelka ~{strelka_vcf}
  >>>

  output {
    File combined_vcf = outfile
    File combined_vcf_tbi = outfile + ".tbi"
  }
}
