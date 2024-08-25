version 1.0

task splitIntervalListToBed {
  input {
    File interval_list
    Int scatter_count
  }

  Int space_needed_gb = 10 + round(size(interval_list, "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "6GB"
    docker: "mgibio/cle:v1.4.2"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /usr/bin/perl /usr/bin/split_interval_list_to_bed_helper.pl OUTPUT="$PWD" INPUT=~{interval_list} SCATTER_COUNT=~{scatter_count}
  >>>

  output {
    Array[File] split_beds = glob("*.interval.bed")
  }
}

workflow wf {
  input {
    File interval_list
    Int scatter_count
  }

  call splitIntervalListToBed {
    input:
    interval_list=interval_list,
    scatter_count=scatter_count
  }
}
