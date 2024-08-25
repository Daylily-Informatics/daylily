version 1.0

task duphold {
  input {
    File bam
    String output_vcf_name = "duphold_annotated.vcf"
    File reference
    File reference_fai
    File reference_dict
    File? snps_vcf
    File sv_vcf
  }

  Int cores = 2
  Int space_needed_gb = 10
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "10GB"
    docker: "mgibio/duphold-cwl:0.1.5"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /usr/bin/duphold --threads ~{cores} \
    --bam ~{bam} \
    --output ~{output_vcf_name} \
    --fasta ~{reference} \
    ~{if defined(snps_vcf) then "--snp " + snps_vcf else ""} \
    --vcf ~{sv_vcf}
  >>>

  output {
    File annotated_sv_vcf = output_vcf_name
  }
}
