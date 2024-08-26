version 1.0

import "../types.wdl"
import "../subworkflows/generate_fda_metrics_for_bam_or_fastqs.wdl" as generate_fda_metrics_for_bam_or_fastqs
import "../tools/cram_to_bam.wdl" as cram_to_bam

workflow generateFdaMetrics {
  input {
    File reference
    File reference_index
    File reference_dict

    Array[SequenceData] unaligned_normal_dna
    Array[SequenceData] unaligned_tumor_dna
    Array[SequenceData] unaligned_tumor_rna
    File aligned_normal_dna
    File aligned_normal_dna_index
    File aligned_tumor_dna
    File aligned_tumor_dna_index
    File aligned_tumor_rna

    QCMetrics normal_qc_metrics
    File normal_duplication_metrics
    QCMetrics tumor_qc_metrics
    File tumor_duplication_metrics
    File rna_metrics

    String? reference_genome
    String? dna_sequencing_platform
    String? dna_sequencing_instrument
    String? dna_sequencing_kit
    String? dna_sequencing_type
    String? dna_single_or_paired_end
    String? normal_dna_spike_in_error_rate
    String? tumor_dna_spike_in_error_rate
    String? normal_dna_total_dna
    String? tumor_dna_total_dna
    String? normal_dna_sample_name
    String? tumor_dna_sample_name
    String? rna_sequencing_platform
    String? rna_sequencing_instrument
    String? rna_sequencing_kit
    String? rna_sequencing_type
    String? rna_single_or_paired_end
    String? rna_spike_in_error_rate
    String? rna_total_rna
    String? rna_rin_score
    String? rna_freq_normalization_method
    String? rna_annotation_file
    String? rna_sample_name
  }

  call generate_fda_metrics_for_bam_or_fastqs.generateFdaMetricsForBamOrFastqs as unaligned_normal_dna_fda_metrics {
    input:
      unaligned_sequence = unaligned_normal_dna,
      aligned_data = [],
      reference = reference,
      reference_index = reference_index,
      stats_output_name = "normal_dna",
      md5sum_output_name = "normal_dna_unaligned",
      table_file_name = "unaligned_normal_dna_table1.csv",
      table_num = "table1",
      is_aligned = false,
      sample_name = normal_dna_sample_name,
      sequencing_platform = dna_sequencing_platform,
      sequencing_instrument = dna_sequencing_instrument,
      sequencing_kit = dna_sequencing_kit,
      single_or_paired_end = dna_single_or_paired_end,
      sequencing_type = dna_sequencing_type,
      spike_in_error_rate = normal_dna_spike_in_error_rate
  }

  call generate_fda_metrics_for_bam_or_fastqs.generateFdaMetricsForBamOrFastqs as unaligned_tumor_dna_fda_metrics {
    input:
      unaligned_sequence = unaligned_tumor_dna,
      aligned_data = [],
      reference = reference,
      reference_index = reference_index,
      stats_output_name = "tumor_dna",
      md5sum_output_name = "tumor_dna_unaligned",
      table_file_name = "unaligned_tumor_dna_table1.csv",
      table_num = "table1",
      is_aligned = false,
      sample_name = tumor_dna_sample_name,
      sequencing_platform = dna_sequencing_platform,
      sequencing_instrument = dna_sequencing_instrument,
      sequencing_kit = dna_sequencing_kit,
      single_or_paired_end = dna_single_or_paired_end,
      sequencing_type = dna_sequencing_type,
      spike_in_error_rate = tumor_dna_spike_in_error_rate
  }

  call generate_fda_metrics_for_bam_or_fastqs.generateFdaMetricsForBamOrFastqs as unaligned_tumor_rna_fda_metrics {
    input:
      unaligned_sequence = unaligned_tumor_rna,
      aligned_data = [],
      reference = reference,
      reference_index = reference_index,
      stats_output_name = "tumor_rna",
      md5sum_output_name = "tumor_rna_unaligned",
      table_file_name = "unaligned_tumor_rna_table1.csv",
      table_num = "table1",
      is_aligned = false,
      sample_name = rna_sample_name,
      sequencing_platform = rna_sequencing_platform,
      sequencing_instrument = rna_sequencing_instrument,
      sequencing_kit = rna_sequencing_kit,
      single_or_paired_end = rna_single_or_paired_end,
      sequencing_type = rna_sequencing_type,
      spike_in_error_rate = rna_spike_in_error_rate
  }

  call cram_to_bam.cramToBam as aligned_normal_dna_cram_to_bam {
    input:
      reference = reference,
      reference_index = reference_index,
      reference_dict = reference_dict,
      cram = aligned_normal_dna,
      cram_index = aligned_normal_dna_index
  }

  call generate_fda_metrics_for_bam_or_fastqs.generateFdaMetricsForBamOrFastqs as aligned_normal_dna_fda_metrics {
    input:
      unaligned_sequence = [],
      aligned_data = [aligned_normal_dna_cram_to_bam.bam],
      reference = reference,
      reference_index = reference_index,
      stats_output_name = "normal_dna",
      md5sum_output_name = "normal_dna_aligned",
      table_file_name = "aligned_normal_dna_table2.csv",
      table_num = "table2",
      is_aligned = true,
      alignment_summary_metrics = normal_qc_metrics.alignment_summary_metrics,
      duplication_metrics = normal_duplication_metrics,
      insert_size_metrics = normal_qc_metrics.insert_size_metrics,
      hs_metrics = normal_qc_metrics.hs_metrics,
      flagstat = normal_qc_metrics.flagstats,
      sample_name = normal_dna_sample_name,
      sequencing_platform = dna_sequencing_platform,
      sequencing_instrument = dna_sequencing_instrument,
      sequencing_kit = dna_sequencing_kit,
      single_or_paired_end = dna_single_or_paired_end,
      spike_in_error_rate = normal_dna_spike_in_error_rate,
      table_source = "Normal sample",
      total_dna = normal_dna_total_dna,
      reference_genome = reference_genome
  }

  call cram_to_bam.cramToBam as aligned_tumor_dna_cram_to_bam {
    input:
      reference = reference,
      reference_index = reference_index,
      reference_dict = reference_dict,
      cram = aligned_tumor_dna,
      cram_index = aligned_tumor_dna_index
  }

  call generate_fda_metrics_for_bam_or_fastqs.generateFdaMetricsForBamOrFastqs as aligned_tumor_dna_fda_metrics {
    input:
      unaligned_sequence = [],
      aligned_data = [aligned_tumor_dna_cram_to_bam.bam],
      reference = reference,
      reference_index = reference_index,
      stats_output_name = "tumor_dna",
      md5sum_output_name = "tumor_dna_aligned",
      table_file_name = "aligned_tumor_dna_table2.csv",
      table_num = "table2",
      is_aligned = true,
      alignment_summary_metrics = tumor_qc_metrics.alignment_summary_metrics,
      duplication_metrics = tumor_duplication_metrics,
      insert_size_metrics = tumor_qc_metrics.insert_size_metrics,
      hs_metrics = tumor_qc_metrics.hs_metrics,
      flagstat = tumor_qc_metrics.flagstats,
      sample_name = tumor_dna_sample_name,
      sequencing_platform = dna_sequencing_platform,
      sequencing_instrument = dna_sequencing_instrument,
      sequencing_kit = dna_sequencing_kit,
      single_or_paired_end = dna_single_or_paired_end,
      spike_in_error_rate = tumor_dna_spike_in_error_rate,
      table_source = "Tumor sample",
      total_dna = tumor_dna_total_dna,
      reference_genome = reference_genome
  }

  call generate_fda_metrics_for_bam_or_fastqs.generateFdaMetricsForBamOrFastqs as aligned_tumor_rna_fda_metrics {
    input:
      unaligned_sequence = [],
      aligned_data = [aligned_tumor_rna],
      reference = reference,
      reference_index = reference_index,
      stats_output_name = "tumor_rna",
      md5sum_output_name = "tumor_rna_aligned",
      table_file_name = "aligned_tumor_rna_table3.csv",
      table_num = "table3",
      is_aligned = true,
      rna_metrics = rna_metrics,
      unaligned_rna_table = unaligned_tumor_rna_fda_metrics.metrics.table,
      sample_name = rna_sample_name,
      sequencing_platform = rna_sequencing_platform,
      sequencing_instrument = rna_sequencing_instrument,
      sequencing_kit = rna_sequencing_kit,
      single_or_paired_end = rna_single_or_paired_end,
      spike_in_error_rate = rna_spike_in_error_rate,
      table_source = "Tumor sample",
      total_rna = rna_total_rna,
      rin_score = rna_rin_score,
      freq_normalization_method = rna_freq_normalization_method,
      annotation_file = rna_annotation_file,
      reference_genome = reference_genome
  }

  output {
    FdaMetrics unaligned_normal_dna_metrics = unaligned_normal_dna_fda_metrics.metrics
    FdaMetrics unaligned_tumor_dna_metrics = unaligned_tumor_dna_fda_metrics.metrics
    FdaMetrics unaligned_tumor_rna_metrics = unaligned_tumor_rna_fda_metrics.metrics
    FdaMetrics aligned_normal_dna_metrics = aligned_normal_dna_fda_metrics.metrics
    FdaMetrics aligned_tumor_dna_metrics = aligned_tumor_dna_fda_metrics.metrics
    FdaMetrics aligned_tumor_rna_metrics = aligned_tumor_rna_fda_metrics.metrics
  }
}
