version 1.0

task filterVcfSomaticLlr {
  input {
    File vcf
    String tumor_sample_name
    String normal_sample_name
    Float threshold = 5
    Float tumor_purity = 1
    Float normal_contamination_rate = 0
  }

  Int space_needed_gb = 10 + round(size(vcf, "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/somatic-llr-filter:v0.4.3"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "somatic_llr_filtered.vcf"
  command <<<
    /opt/conda/bin/python3 /usr/bin/somatic_llr_filter.py \
    --normal-contamination-rate ~{normal_contamination_rate} \
    --tumor-purity ~{tumor_purity} \
    --tumor-sample-name ~{tumor_sample_name} \
    --llr-threshold ~{threshold} \
    --overwrite \
    --normal-sample-name ~{normal_sample_name} \
    ~{vcf} ~{outfile} \
  >>>

  output {
    File somatic_llr_filtered_vcf = outfile
  }
}
