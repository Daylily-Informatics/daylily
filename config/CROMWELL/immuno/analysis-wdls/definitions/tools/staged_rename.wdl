version 1.0

task stagedRename {
  input {
    File original
    String name
  }

  Int space_needed_gb = 10 + round(size(original, "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    cpu: 1
    docker: "ubuntu:bionic"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /bin/mv ~{original} ~{name}
  >>>

  output {
    File replacement = name
  }
}
