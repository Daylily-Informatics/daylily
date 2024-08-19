version 1.0

workflow AlignFastq {
  input {
    File fastq1
    File fastq2
    File reference_fasta
    File reference_fasta_index
    File reference_fasta_dict
    String docker_image
    String partition
    Int cpus
    String memory
  }

  call AlignFastq {
    input:
      fastq1 = fastq1,
      fastq2 = fastq2,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
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
    File reference_fasta_index
    File reference_fasta_dict
    String docker_image
    String partition
    Int cpus
    String memory
  }

  command {
    docker run ${docker_image} \
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
  }
}