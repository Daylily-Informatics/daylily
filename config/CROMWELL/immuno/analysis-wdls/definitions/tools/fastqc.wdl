version 1.0

task fastqc {
  input {
    Array[File] files
  }


  Int space_needed_gb = 10 + round(size(files, "GB"))
  runtime {
    memory: "32GB"
    docker: "mgibio/fastqc:0.11.9"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    mkdir outputs
    /usr/bin/perl /usr/local/bin/FastQC/fastqc -q -o outputs ~{sep=" " files}
  >>>

  output {
    Array[File] fastqc_data = glob("outputs/*.zip")
  }
}
