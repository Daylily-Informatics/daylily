version 1.0

task collectWgsMetrics {
  input {
    String sample_name = "final"
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
    File? intervals
    Int? minimum_mapping_quality
    Int? minimum_base_quality
  }

  Float bam_size = size([bam, bam_bai], "GB")
  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Float intervals_size = size(intervals, "GB")
  Int space_needed_gb = 10 + round(bam_size + reference_size + intervals_size)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "48GB"
    docker: "broadinstitute/picard:2.23.6"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outname = sample_name + ".WgsMetrics.txt"
  command <<<
    /usr/bin/java -Xmx32g -jar /usr/picard/picard.jar CollectWgsMetrics \
    O=~{outname} \
    I=~{bam} R=~{reference} \
    ~{if defined(intervals) then "INTERVALS=" + select_first([intervals]) else ""} \
    ~{if defined(minimum_mapping_quality) then "MINIMUM_MAPPING_QUALITY=" + select_first([minimum_mapping_quality]) else ""} \
    ~{if defined(minimum_base_quality) then "MINIMUM_BASE_QUALITY=" + select_first([minimum_base_quality]) else ""}
  >>>

  output {
    File wgs_metrics = outname
  }
}
