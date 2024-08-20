version 1.0

workflow AlignFastq {
  input {
    File fastq1
    File fastq2
    File reference_fasta
    File reference_fasta_dict
    File reference_fasta_amb
    File reference_fasta_ann
    File reference_fasta_bwt
    File reference_fasta_pac
    File reference_fasta_sa
    String docker_image
    String partition
    Int cpu
    Int memory

  }
  
  call AlignFastq {
    input:
      fastq1 = fastq1,
      fastq2 = fastq2,
      reference_fasta = reference_fasta,
      reference_fasta_amb = reference_fasta_amb,
      reference_fasta_ann = reference_fasta_ann,
      reference_fasta_bwt = reference_fasta_bwt,
      reference_fasta_pac = reference_fasta_pac,
      reference_fasta_sa = reference_fasta_sa,
      reference_fasta_dict = reference_fasta_dict,
      docker_image = docker_image,
      partition = partition,
      cpu = cpu,
      memory = memory,
  }

  output {
    File bam = output.bam
  }
}

task AlignFastq {
  input {
    File fastq1
    File fastq2
    File reference_fasta
    File reference_fasta_amb
    File reference_fasta_ann
    File reference_fasta_bwt
    File reference_fasta_pac
    File reference_fasta_sa
    File reference_fasta_dict
    String docker_image
    String partition
    Int cpu
    Int memory
    String project
    String all_partitions
  }

  command {
      bwa mem -t ${cpu} ${reference_fasta} ${fastq1} ${fastq2} | \
      samtools view -Sb - > output.bam
  }

  output {
    File bam = "output.bam"
  }

  runtime {
    cpu: cpu
    memory: memory
    docker: docker_image
    partition: partition
    project: project
    all_partitions: all_partitions
  }
}