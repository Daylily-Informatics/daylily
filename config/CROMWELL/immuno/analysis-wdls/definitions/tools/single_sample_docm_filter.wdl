version 1.0

task singleSampleDocmFilter {
  input {
    File docm_out
  }

  Int space_needed_gb = 10 + round(size(docm_out, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    docker: "mgibio/perl_helper-cwl:1.0.0"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /usr/bin/perl /usr/bin/single_sample_docm_filter.pl ~{docm_out} "$PWD"
  >>>

  output {
    File docm_filter_out = "docm_filter_out.vcf"
  }
}
