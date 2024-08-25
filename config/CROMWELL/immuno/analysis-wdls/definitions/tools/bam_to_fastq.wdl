version 1.0

task bamToFastq {
  input { File bam }

  # ran into issue at 3*, bump to 10*
  Int space_needed_gb = 10 + round(10*size(bam, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/rnaseq:1.0.0"
    cpu: 1
    memory: "6GB"
    bootDiskSizeGb: 25
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar SamToFastq VALIDATION_STRINGENCY=SILENT \
        F=read1.fastq F2=read2.fastq I=~{bam}
  >>>

  output {
    File fastq1 = "read1.fastq"
    File fastq2 = "read2.fastq"
    Array[File] fastqs = ["read1.fastq", "read2.fastq"]
  }
}

workflow wf {
  input { File bam }
  call bamToFastq { input: bam=bam }
}
