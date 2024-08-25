version 1.0

task catAll {
  input {
    Array[File] region_pindel_outs
  }

  Int space_needed_gb = 10 + round(size(region_pindel_outs, "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    docker: "ubuntu:xenial"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /bin/grep "ChrID" ~{sep=" " region_pindel_outs} > all_region_pindel.head
  >>>

  output {
    File all_region_pindel_head = "all_region_pindel.head"
  }
}
