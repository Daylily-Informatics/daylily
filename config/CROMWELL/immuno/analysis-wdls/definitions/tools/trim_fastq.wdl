version 1.0

task trimFastq {
  input {
    File adapters
    String adapter_trim_end
    Int adapter_min_overlap
    Int max_uncalled
    Int min_readlength
    File reads1
    File reads2
  }

  Int cores = 4
  Int space_needed_gb = 10 + ceil(size(adapters, "GB") + 10*size([reads1, reads2], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "16GB"
    bootDiskSizeGb: 25
    cpu: 4
    docker: "mgibio/bisulfite:v1.4"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /opt/flexbar/flexbar --target trimmed_read --threads ~{cores} \
        --adapters ~{adapters} \
        --adapter-trim-end ~{adapter_trim_end} \
        --adapter-min-overlap ~{adapter_min_overlap} \
        --max-uncalled ~{max_uncalled} \
        --min-read-length ~{min_readlength} \
        --reads ~{reads1} --reads2 ~{reads2}
  >>>

  output {
    File fastq1 = "trimmed_read_1.fastq"
    File fastq2 = "trimmed_read_2.fastq"
    Array[File] fastqs = glob("trimmed_read_*.fastq")
  }
}

workflow wf {
  input {
    File adapters
    String adapter_trim_end
    Int adapter_min_overlap
    Int max_uncalled
    Int min_readlength
    File reads1
    File reads2
  }
  call trimFastq {
    input:
    adapters=adapters,
    adapter_trim_end=adapter_trim_end,
    adapter_min_overlap=adapter_min_overlap,
    max_uncalled=max_uncalled,
    min_readlength=min_readlength,
    reads1=reads1,
    reads2=reads2
  }
}
