version 1.0

task pvacfuse {
  input {
    File input_fusions_zip
    File? peptide_fasta
    String sample_name
    Array[String] alleles
    Array[String] prediction_algorithms
    Array[Int]? epitope_lengths_class_i
    Array[Int]? epitope_lengths_class_ii
    Int? binding_threshold
    Int? percentile_threshold
    Int? iedb_retries
    Boolean keep_tmp_files = false
    String? net_chop_method  # enum [cterm 20s]
    Boolean netmhc_stab = false
    String? top_score_metric  # enum [lowest, median]
    Float? net_chop_threshold
    Boolean run_reference_proteome_similarity = false
    String? additional_report_columns  # enum [sample_name]
    Int? fasta_size
    Int? downstream_sequence_length
    Boolean exclude_nas = false
    Int n_threads = 8
    Boolean allele_specific_binding_thresholds = false
    Int? aggregate_inclusion_binding_threshold
    Array[String]? problematic_amino_acids
    File? star_fusion_file
    Int? read_support
    Float? expn_val
  }

  Int space_needed_gb = 10 + round(size([input_fusions_zip], "GB") * 3)
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "griffithlab/pvactools:4.2.0"
    memory: "16GB"
    cpu: n_threads
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  # explicit typing required, don't inline
  Array[Int] epitope_i = select_first([epitope_lengths_class_i, []])
  Array[Int] epitope_ii = select_first([epitope_lengths_class_ii, []])
  Array[String] problematic_aa = select_first([problematic_amino_acids, []])
  command <<<
    mkdir agfusion_dir && unzip -qq ~{input_fusions_zip} -d agfusion_dir

    ln -s "$TMPDIR" /tmp/pvacfuse && export TMPDIR=/tmp/pvacfuse && \
    /usr/local/bin/pvacfuse run --iedb-install-directory /opt/iedb \
    agfusion_dir ~{sample_name} \
    ~{sep="," alleles} \
    ~{sep=" " prediction_algorithms} \
    pvacfuse_predictions \
    ~{if length(epitope_i) > 0 then "-e1 " else ""} ~{sep="," epitope_i} \
    ~{if length(epitope_ii) > 0 then "-e2 " else ""} ~{sep="," epitope_ii} \
    ~{if defined(binding_threshold) then "-b ~{binding_threshold}" else ""} \
    ~{if defined(percentile_threshold) then "--percentile-threshold ~{percentile_threshold}" else ""} \
    ~{if allele_specific_binding_thresholds then "--allele-specific-binding-thresholds" else ""} \
    ~{if defined(aggregate_inclusion_binding_threshold) then "--aggregate-inclusion-binding-threshold ~{aggregate_inclusion_binding_threshold}" else ""} \
    ~{if defined(iedb_retries) then "-r ~{iedb_retries}" else ""} \
    ~{if keep_tmp_files then "-k" else ""} \
    ~{if defined(net_chop_method) then "--net-chop-method ~{net_chop_method}" else ""} \
    ~{if netmhc_stab then "--netmhc-stab" else ""} \
    ~{if defined(top_score_metric) then "-m ~{top_score_metric}" else ""} \
    ~{if defined(top_score_metric) then "-m ~{top_score_metric}" else ""} \
    ~{if defined(net_chop_threshold) then "--net-chop-threshold ~{net_chop_threshold}" else ""} \
    ~{if run_reference_proteome_similarity then "--run-reference-proteome-similarity" else ""} \
    ~{if defined(peptide_fasta) then "--peptide-fasta ~{peptide_fasta}" else ""} \
    ~{if defined(additional_report_columns) then "-m ~{additional_report_columns}" else ""} \
    ~{if defined(fasta_size) then "-s ~{fasta_size}" else ""} \
    ~{if defined(downstream_sequence_length) then "-d ~{downstream_sequence_length}" else ""} \
    ~{if exclude_nas then "--exclude-NAs" else ""} \
    ~{if length(problematic_aa) > 0 then "--problematic-amino-acids" else ""} ~{sep="," problematic_aa} \
    ~{if defined(star_fusion_file) then "--starfusion-file ~{star_fusion_file}" else ""} \
    ~{if defined(read_support) then "--read-support ~{read_support}" else ""} \
    ~{if defined(expn_val) then "--expn-val ~{expn_val}" else ""} \
    --n-threads ~{n_threads}
  >>>

  output {
    File? mhc_i_all_epitopes = "pvacfuse_predictions/MHC_Class_I/~{sample_name}.all_epitopes.tsv"
    File? mhc_i_aggregated_report = "pvacfuse_predictions/MHC_Class_I/~{sample_name}.all_epitopes.aggregated.tsv"
    File? mhc_i_filtered_epitopes = "pvacfuse_predictions/MHC_Class_I/~{sample_name}.filtered.tsv"
    File? mhc_ii_all_epitopes = "pvacfuse_predictions/MHC_Class_II/~{sample_name}.all_epitopes.tsv"
    File? mhc_ii_aggregated_report = "pvacfuse_predictions/MHC_Class_II/~{sample_name}.all_epitopes.aggregated.tsv"
    File? mhc_ii_filtered_epitopes = "pvacfuse_predictions/MHC_Class_II/~{sample_name}.filtered.tsv"
    File? combined_all_epitopes = "pvacfuse_predictions/combined/~{sample_name}.all_epitopes.tsv"
    File? combined_aggregated_report = "pvacfuse_predictions/combined/~{sample_name}.all_epitopes.aggregated.tsv"
    File? combined_filtered_epitopes = "pvacfuse_predictions/combined/~{sample_name}.filtered.tsv"

    # glob documentation
    # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md#globs
    Array[File] mhc_i = glob("pvacfuse_predictions/MHC_Class_I/*")
    Array[File] mhc_ii = glob("pvacfuse_predictions/MHC_Class_II/*") 
    Array[File] combined = glob("pvacfuse_predictions/combined/*")

  }
}
