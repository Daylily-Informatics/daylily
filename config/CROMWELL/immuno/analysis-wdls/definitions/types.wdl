version 1.0

#
# Input Types
#

struct Sequence {
  File? bam
  File? fastq1
  File? fastq2
}
# assume either bam or fastqs defined
struct SequenceData {
  Sequence sequence
  String? readgroup
}

struct TrimmingOptions {
  File adapters
  Int min_overlap
}

struct LabelledFile {
  File file
  String label
}

# ---- vep_custom_annotation ----
struct Info {
  File file
  Array[File]? secondary_files
  String data_format  # enum, ['bed', 'gff', 'gtf', 'vcf', 'bigwig']
  String name
  Array[String]? vcf_fields
  Boolean? gnomad_filter
  Boolean check_existing
}

struct VepCustomAnnotation {
  Boolean force_report_coordinates
  String method  # enum, ['exact', 'overlap']
  Info annotation
}

#
# Output Types
#

struct QCMetrics {
  File insert_size_metrics
  File insert_size_histogram
  File alignment_summary_metrics
  File? hs_metrics
  File? gc_bias_metrics
  File? gc_bias_metrics_chart
  File? gc_bias_metrics_summary
  File? wgs_metrics
  Array[File] per_target_coverage_metrics
  Array[File] per_target_hs_metrics
  Array[File] per_base_coverage_metrics
  Array[File] per_base_hs_metrics
  Array[File] summary_hs_metrics
  File flagstats
  File? verify_bam_id_metrics
  File? verify_bam_id_depth
  File? bamcoverage_bigwig
}
