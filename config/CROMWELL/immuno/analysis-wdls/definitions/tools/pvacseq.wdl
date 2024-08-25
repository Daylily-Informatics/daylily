version 1.0

task pvacseq {
  input {
    Int n_threads = 8
    File input_vcf
    File input_vcf_tbi
    String sample_name
    Array[String] alleles
    Array[String] prediction_algorithms
    File? peptide_fasta

    Array[Int]? epitope_lengths_class_i
    Array[Int]? epitope_lengths_class_ii
    Int? binding_threshold
    Int? percentile_threshold
    Int? iedb_retries

    String? normal_sample_name
    String? net_chop_method  # enum [cterm , 20s]
    String? top_score_metric  # enum [lowest, median]
    Float? net_chop_threshold
    String? additional_report_columns  # enum [sample_name]
    Int? fasta_size
    Int? downstream_sequence_length
    Boolean exclude_nas = false
    File? phased_proximal_variants_vcf
    File? phased_proximal_variants_vcf_tbi
    Float? minimum_fold_change
    Int? normal_cov
    Int? tdna_cov
    Int? trna_cov
    Float? normal_vaf
    Float? tdna_vaf
    Float? trna_vaf
    Float? expn_val
    Int? maximum_transcript_support_level  # enum [1, 2, 3, 4, 5]
    Int? aggregate_inclusion_binding_threshold
    Array[String]? problematic_amino_acids
    Float? anchor_contribution_threshold

    Boolean allele_specific_binding_thresholds = false
    Boolean keep_tmp_files = false
    Boolean netmhc_stab = false
    Boolean run_reference_proteome_similarity = false
    Boolean allele_specific_anchors = false

    Float? tumor_purity
  }

  Float input_size = size([input_vcf, input_vcf_tbi], "GB")
  Float phased_variants_size = size([phased_proximal_variants_vcf, phased_proximal_variants_vcf_tbi], "GB")
  Int space_needed_gb = 10 + round(input_size + phased_variants_size)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "16GB"
    cpu: n_threads
    docker: "griffithlab/pvactools:4.2.0"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  # explicit typing required, don't inline
  Array[Int] epitope_i = select_first([epitope_lengths_class_i, []])
  Array[Int] epitope_ii = select_first([epitope_lengths_class_ii, []])
  Array[String] problematic_aa = select_first([problematic_amino_acids, []])
  command <<<
    # touch each tbi to ensure they have a timestamp after the vcf
    touch ~{phased_proximal_variants_vcf_tbi}
    touch ~{input_vcf_tbi}

    ln -s "$TMPDIR" /tmp/pvacseq && export TMPDIR=/tmp/pvacseq && \
    /usr/local/bin/pvacseq run --iedb-install-directory /opt/iedb \
    --pass-only \
    ~{if defined(tumor_purity) then "--tumor-purity " + select_first([tumor_purity]) else ""} \
    ~{if length(epitope_i ) > 0 then "-e1 " else ""} ~{sep="," epitope_i} \
    ~{if length(epitope_ii) > 0 then "-e2 " else ""} ~{sep="," epitope_ii} \
    ~{if defined(binding_threshold) then "-b ~{binding_threshold}" else ""} \
    ~{if defined(percentile_threshold) then "--percentile-threshold ~{percentile_threshold}" else ""} \
    ~{if allele_specific_binding_thresholds then "--allele-specific-binding-thresholds" else ""} \
    ~{if defined(aggregate_inclusion_binding_threshold) then "--aggregate-inclusion-binding-threshold ~{aggregate_inclusion_binding_threshold}" else ""} \
    ~{if defined(iedb_retries) then "-r ~{iedb_retries}" else ""} \
    ~{if keep_tmp_files then "-k" else ""} \
    ~{if defined(normal_sample_name) then "--normal-sample-name ~{normal_sample_name}" else ""} \
    ~{if defined(net_chop_method) then "--net-chop-method ~{net_chop_method}" else ""} \
    ~{if netmhc_stab then "--netmhc-stab" else ""} \
    ~{if run_reference_proteome_similarity then "--run-reference-proteome-similarity" else ""} \
    ~{if defined(peptide_fasta) then "--peptide-fasta ~{peptide_fasta}" else ""} \
    ~{if defined(top_score_metric) then "-m ~{top_score_metric}" else ""} \
    ~{if defined(net_chop_threshold) then "--net-chop-threshold ~{net_chop_threshold}" else ""} \
    ~{if defined(additional_report_columns) then "-m ~{additional_report_columns}" else ""} \
    ~{if defined(fasta_size) then "-s ~{fasta_size}" else ""} \
    ~{if defined(downstream_sequence_length) then "-d ~{downstream_sequence_length}" else ""} \
    ~{if exclude_nas then "--exclude-NAs" else ""} \
    ~{if defined(phased_proximal_variants_vcf) then "-p ~{phased_proximal_variants_vcf}" else ""} \
    ~{if defined(minimum_fold_change) then "-c ~{minimum_fold_change}" else ""} \
    ~{if defined(normal_cov) then "--normal-cov ~{normal_cov}" else ""} \
    ~{if defined(tdna_cov) then "--tdna-cov ~{tdna_cov}" else ""} \
    ~{if defined(trna_cov) then "--trna-cov ~{trna_cov}" else ""} \
    ~{if defined(normal_vaf) then "--normal-vaf ~{normal_vaf}" else ""} \
    ~{if defined(tdna_vaf) then "--tdna-vaf ~{tdna_vaf}" else ""} \
    ~{if defined(trna_vaf) then "--trna-vaf ~{trna_vaf}" else ""} \
    ~{if defined(expn_val) then "--expn-val ~{expn_val}" else ""} \
    ~{if defined(maximum_transcript_support_level) then "--maximum-transcript-support-level ~{maximum_transcript_support_level}" else ""} \
    ~{if length(problematic_aa) > 0 then "--problematic-amino-acids" else ""} ~{sep="," problematic_aa} \
    ~{if allele_specific_anchors then "--allele-specific-anchors" else ""} \
    ~{if defined(anchor_contribution_threshold) then "--anchor-contribution-threshold ~{anchor_contribution_threshold}" else ""} \
    --n-threads ~{n_threads} \
    ~{input_vcf} ~{sample_name} ~{sep="," alleles} ~{sep=" " prediction_algorithms} \
    pvacseq_predictions
  >>>

  output {
    File? mhc_i_all_epitopes = "pvacseq_predictions/MHC_Class_I/~{sample_name}.all_epitopes.tsv"
    File? mhc_i_aggregated_report = "pvacseq_predictions/MHC_Class_I/~{sample_name}.all_epitopes.aggregated.tsv"
    File? mhc_i_filtered_epitopes = "pvacseq_predictions/MHC_Class_I/~{sample_name}.filtered.tsv"
    File? mhc_i_aggregated_metrics_file = "pvacseq_predictions/MHC_Class_I/" + sample_name + ".all_epitopes.aggregated.metrics.json"
    File? mhc_ii_all_epitopes = "pvacseq_predictions/MHC_Class_II/~{sample_name}.all_epitopes.tsv"
    File? mhc_ii_aggregated_report = "pvacseq_predictions/MHC_Class_II/~{sample_name}.all_epitopes.aggregated.tsv"
    File? mhc_ii_filtered_epitopes = "pvacseq_predictions/MHC_Class_II/~{sample_name}.filtered.tsv"
    File? mhc_ii_aggregated_metrics_file = "pvacseq_predictions/MHC_Class_II/" + sample_name + ".all_epitopes.aggregated.metrics.json"
    File? combined_all_epitopes = "pvacseq_predictions/combined/~{sample_name}.all_epitopes.tsv"
    File? combined_aggregated_report = "pvacseq_predictions/combined/~{sample_name}.all_epitopes.aggregated.tsv"
    File? combined_filtered_epitopes = "pvacseq_predictions/combined/~{sample_name}.filtered.tsv"
    File? combined_aggregated_metrics_file = "pvacseq_predictions/combined/" + sample_name + ".all_epitopes.aggregated.metrics.json"
    
    # glob documentations
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#globs
    Array[File] mhc_i = glob("pvacseq_predictions/MHC_Class_I/*")
    Array[File] mhc_ii = glob("pvacseq_predictions/MHC_Class_II/*")
    Array[File] combined = glob("pvacseq_predictions/combined/*")

  }
}

workflow wf {
  input {
    Int? n_threads
    File input_vcf
    File input_vcf_tbi
    String sample_name
    Array[String] alleles
    Array[String] prediction_algorithms

    Array[Int]? epitope_lengths_class_i
    Array[Int]? epitope_lengths_class_ii
    Int? binding_threshold
    Int? percentile_threshold
    Int? iedb_retries

    String? normal_sample_name
    String? net_chop_method  # enum [cterm , 20s]
    String? top_score_metric  # enum [lowest, median]
    Float? net_chop_threshold
    String? additional_report_columns  # enum [sample_name]
    Int? fasta_size
    Int? downstream_sequence_length
    Boolean exclude_nas = false
    File? phased_proximal_variants_vcf
    File? phased_proximal_variants_vcf_tbi
    Float? minimum_fold_change
    Int? normal_cov
    Int? tdna_cov
    Int? trna_cov
    Float? normal_vaf
    Float? tdna_vaf
    Float? trna_vaf
    Float? expn_val
    String? maximum_transcript_support_level  # enum [1, 2, 3, 4, 5]
    Int? aggregate_inclusion_binding_threshold
    Array[String]? problematic_amino_acids
    Float? anchor_contribution_threshold

    Boolean? allele_specific_binding_thresholds
    Boolean? keep_tmp_files
    Boolean? netmhc_stab
    Boolean? run_reference_proteome_similarity
    Boolean allele_specific_anchors = false
  }
  call pvacseq {
    input:
    n_threads=n_threads,
    input_vcf=input_vcf,
    input_vcf_tbi=input_vcf_tbi,
    sample_name=sample_name,
    alleles=alleles,
    prediction_algorithms=prediction_algorithms,
    epitope_lengths_class_i=epitope_lengths_class_i,
    epitope_lengths_class_ii=epitope_lengths_class_ii,
    binding_threshold=binding_threshold,
    percentile_threshold=percentile_threshold,
    aggregate_inclusion_binding_threshold=aggregate_inclusion_binding_threshold,
    iedb_retries=iedb_retries,
    normal_sample_name=normal_sample_name,
    net_chop_method=net_chop_method,
    top_score_metric=top_score_metric,
    net_chop_threshold=net_chop_threshold,
    additional_report_columns=additional_report_columns,
    fasta_size=fasta_size,
    downstream_sequence_length=downstream_sequence_length,
    exclude_nas=exclude_nas,
    phased_proximal_variants_vcf=phased_proximal_variants_vcf,
    phased_proximal_variants_vcf_tbi=phased_proximal_variants_vcf_tbi,
    minimum_fold_change=minimum_fold_change,
    normal_cov=normal_cov,
    tdna_cov=tdna_cov,
    trna_cov=trna_cov,
    normal_vaf=normal_vaf,
    tdna_vaf=tdna_vaf,
    trna_vaf=trna_vaf,
    expn_val=expn_val,
    maximum_transcript_support_level=maximum_transcript_support_level,
    allele_specific_binding_thresholds=allele_specific_binding_thresholds,
    problematic_amino_acids=problematic_amino_acids,
    allele_specific_anchors=allele_specific_anchors,
    anchor_contribution_threshold=anchor_contribution_threshold,
    keep_tmp_files=keep_tmp_files,
    netmhc_stab=netmhc_stab,
    run_reference_proteome_similarity=run_reference_proteome_similarity
  }
}
