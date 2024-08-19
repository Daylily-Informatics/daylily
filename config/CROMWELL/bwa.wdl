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
