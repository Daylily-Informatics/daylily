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
    Int cpus
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
      cpus = cpus,
      memory = memory
  }

  output {
    File bam = AlignFastq.bam
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
    Int cpus
    Int memory
  }

  command {
      bwa mem -t ${cpus} ${reference_fasta} ${fastq1} ${fastq2} | \
      samtools view -Sb - > output.bam
  }

  output {
    File bam = "output.bam"
  }

  runtime {
    cpu: cpus
    memory: memory
    docker: docker_image
    partition: partition
    entrypoint: "/bin/sh"
  }
}