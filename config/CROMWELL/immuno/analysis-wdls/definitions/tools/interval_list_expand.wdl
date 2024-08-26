version 1.0

task intervalListExpand {
  input {
    File interval_list
    Int roi_padding
  }

  Int space_needed_gb = 10 + round(size(interval_list, "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    docker: "broadinstitute/picard:2.23.6"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String output_file = basename(interval_list, ".interval_list") + ".expanded.interval_list"
  command <<<
    /usr/bin/java -Xmx3g -jar /usr/picard/picard.jar IntervalListTools OUTPUT=~{output_file} UNIQUE=TRUE INPUT=~{interval_list} PADDING=~{roi_padding}
  >>>

  output {
    File expanded_interval_list = output_file
  }
}

workflow wf {
  input {
    File interval_list
    Int roi_padding
  }
  call intervalListExpand {
    input: interval_list=interval_list, roi_padding=roi_padding
  }
}
