version 1.0

task fpFilter {
  input {
    File reference
    File reference_fai
    File reference_dict

    File bam
    File vcf

    String output_vcf_basename = "fpfilter"
    String sample_name = "TUMOR"
    Float fp_min_var_freq = 0.05
  }

  Int space_needed_gb = 10 + round(size(vcf, "GB")*2 + size([reference, reference_fai, reference_dict, bam], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "6GB"
    bootDiskSizeGb: 25
    docker: "mgibio/fp_filter-cwl:1.0.1"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String output_vcf = output_vcf_basename + ".vcf"
  command <<<
    /usr/bin/perl /usr/bin/fpfilter.pl --bam-readcount /usr/bin/bam-readcount --samtools /opt/samtools/bin/samtools --output ~{output_vcf} --reference ~{reference} --bam-file ~{bam} --vcf-file ~{vcf} --sample ~{sample_name} --min-var-freq ~{fp_min_var_freq}
  >>>

  output {
    File filtered_vcf = output_vcf
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict

    File bam
    File vcf

    String sample_name = "TUMOR"
    String output_vcf_basename = "fpfilter"
  }

  call fpFilter {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    vcf=vcf,
    output_vcf_basename=output_vcf_basename,
    sample_name=sample_name
  }
}
