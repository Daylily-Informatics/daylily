version 1.0

task md5sum {
  input {
    Array[File]+ files
    String output_name = basename(files[0])
  }

  Int space_needed_gb = 10 + round(size(files, "GB"))
  runtime {
    memory: "4GB"
    docker: "ubuntu:bionic"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    md5sum ~{sep=" " files} > "~{output_name}.md5"
  >>>

  output {
    File md5_file = "~{output_name}.md5"
  }
}
