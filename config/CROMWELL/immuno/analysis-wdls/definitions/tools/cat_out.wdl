version 1.0

task catOut {
  input {
    Array[File] pindel_outs
  }

  Int space_needed_gb = 10 + round(size(pindel_outs, "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    docker: "ubuntu:xenial"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /bin/cat ~{sep=" " pindel_outs} > "per_chromosome_pindel.out"
  >>>

  output {
    File pindel_out = "per_chromosome_pindel.out"
  }
}
