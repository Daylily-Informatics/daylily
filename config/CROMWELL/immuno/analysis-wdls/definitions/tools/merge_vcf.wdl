version 1.0

task mergeVcf {
  input {
    Array[File] vcfs
    Array[File] vcf_tbis
    String merged_vcf_basename = "merged"
  }

  Int space_needed_gb = 10 + round(2*(size(vcfs, "GB") + size(vcf_tbis, "GB")))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/bcftools-cwl:1.12"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String output_file = merged_vcf_basename + ".vcf.gz"
  command <<<
    /opt/bcftools/bin/bcftools concat --allow-overlaps --remove-duplicates --output-type z -o ~{output_file} ~{sep=" " vcfs}
  >>>

  output {
    File merged_vcf = output_file
  }
}
