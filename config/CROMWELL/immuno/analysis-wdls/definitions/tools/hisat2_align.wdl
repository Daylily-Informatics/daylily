version 1.0

task hisat2Align {
  input {
    File reference_index
    File reference_index_1ht2
    File reference_index_2ht2
    File reference_index_3ht2
    File reference_index_4ht2
    File reference_index_5ht2
    File reference_index_6ht2
    File reference_index_7ht2
    File reference_index_8ht2
    File fastq1
    File fastq2
    String read_group_id
    Array[String] read_group_fields
    String strand = "unstranded"  # [first, second, unstranded]
  }

  Int cores = 16
  Float fastq_size_gb = size([fastq1, fastq2], "GB")
  Float reference_size_gb = size([
    reference_index,
    reference_index_1ht2,
    reference_index_2ht2,
    reference_index_3ht2,
    reference_index_4ht2,
    reference_index_5ht2,
    reference_index_6ht2,
    reference_index_7ht2,
    reference_index_8ht2
  ], "GB")
  Int space_needed_gb = 10 + round(5*fastq_size_gb + reference_size_gb)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "32GB"
    cpu: cores
    bootDiskSizeGb: space_needed_gb
    docker: "mgibio/hisat2-sambamba:0.1"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "aligned.bam"
  Map[String, String] strandness = {
    "first":      "--rna-strandness RF",
    "second":     "--rna-strandness FR",
    "unstranded": ""
  }

  command <<<
    /usr/local/bin/hisat2 \
            ~{strandness[strand]} \
            --rg-id ~{read_group_id} \
            ~{sep=' ' prefix("--rg ", read_group_fields)} \
            -x ~{reference_index} \
            -1 ~{fastq1} \
            -2 ~{fastq2} \
            -p ~{cores} --dta \
        | /usr/local/bin/sambamba view -S -f bam -l 0 /dev/stdin \
        | /usr/local/bin/sambamba sort -t ~{cores} -m 8G -o ~{outfile} /dev/stdin
  >>>

  output {
    File aligned_bam = outfile
  }
}

workflow wf {
  input {
    File reference_index
    File reference_index_1ht2
    File reference_index_2ht2
    File reference_index_3ht2
    File reference_index_4ht2
    File reference_index_5ht2
    File reference_index_6ht2
    File reference_index_7ht2
    File reference_index_8ht2
    File fastq1
    File fastq2
    String read_group_id
    Array[String] read_group_fields
    String? strand
  }
  call hisat2Align {
    input:
    reference_index=reference_index,
    reference_index_1ht2=reference_index_1ht2,
    reference_index_2ht2=reference_index_2ht2,
    reference_index_3ht2=reference_index_3ht2,
    reference_index_4ht2=reference_index_4ht2,
    reference_index_5ht2=reference_index_5ht2,
    reference_index_6ht2=reference_index_6ht2,
    reference_index_7ht2=reference_index_7ht2,
    reference_index_8ht2=reference_index_8ht2,
    fastq1=fastq1,
    fastq2=fastq2,
    read_group_id=read_group_id,
    read_group_fields=read_group_fields,
    strand=strand
  }
}
