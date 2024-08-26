version 1.0

task gatkHaplotypeCaller {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bai
    String emit_reference_confidence  # enum [NONE, BP_RESOLUTION, GVCF]
    Array[String] gvcf_gq_bands
    Array[String] intervals
    File? dbsnp_vcf
    File? dbsnp_vcf_tbi
    File? contamination_fraction  # contains a float. See `tools/freemix.wdl` for context
    Int? max_alternate_alleles
    Int? ploidy
    String? read_filter
    String output_file_name
  }

  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Float bam_size = size([bam, bai], "GB")
  Float vcf_size = size([dbsnp_vcf, dbsnp_vcf_tbi], "GB")
  Int space_needed_gb = 10 + round(reference_size + 2*bam_size + vcf_size)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "18GB"
    docker: "broadinstitute/gatk:4.1.8.1"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  Array[String] pref_bands = prefix("-GQB ", gvcf_gq_bands)
  command <<<
    # requires .bai not .bam.bai
    mv ~{bam} ~{basename(bam)}
    mv ~{bai} ~{basename(basename(bai, ".bai"), ".bam") + ".bai"}

    CONTAMINATION_FILE=~{contamination_fraction}
    if [ -z "$CONTAMINATION_FILE" ]; then
        CONTAMINATION_ARG="--contamination $(head -n 1 $CONTAMINATION_FILE})"
    fi

    # do the task itself
    /gatk/gatk --java-options -Xmx16g HaplotypeCaller \
    -R ~{reference} \
    -I ~{basename(bam)} \
    -ERC ~{emit_reference_confidence} \
    ~{sep=" " pref_bands} \
    -L ~{sep="," intervals} \
    ~{if(defined(dbsnp_vcf)) then "--dbsnp " + dbsnp_vcf else ""} \
    $CONTAMINATION_ARG \
    ~{if(defined(max_alternate_alleles)) then "--max_alternate_alleles " + max_alternate_alleles else ""} \
    ~{if(defined(ploidy)) then "-ploidy " + ploidy else ""} \
    ~{if(defined(read_filter)) then "--read_filter " + read_filter else ""} \
    -O ~{output_file_name}
  >>>

  output {
    File gvcf = output_file_name
    File gvcf_tbi = output_file_name + ".tbi"
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bai
    String emit_reference_confidence  # enum [NONE, BP_RESOLUTION, GVCF]
    Array[String] gvcf_gq_bands
    Array[String] intervals
    File? dbsnp_vcf
    File? dbsnp_vcf_tbi
    File? contamination_fraction
    Int? max_alternate_alleles
    Int? ploidy
    String? read_filter
    String output_file_name
  }
  call gatkHaplotypeCaller {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam,
    bai=bai,
    emit_reference_confidence=emit_reference_confidence,
    gvcf_gq_bands=gvcf_gq_bands,
    intervals=intervals,
    dbsnp_vcf=dbsnp_vcf,
    dbsnp_vcf_tbi=dbsnp_vcf_tbi,
    contamination_fraction=contamination_fraction,
    max_alternate_alleles=max_alternate_alleles,
    ploidy=ploidy,
    read_filter=read_filter,
    output_file_name=output_file_name,
  }
}
