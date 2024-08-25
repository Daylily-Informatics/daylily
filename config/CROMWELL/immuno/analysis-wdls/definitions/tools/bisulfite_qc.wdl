version 1.0

task bisulfiteQc {
  input {
    File vcf
    File bam
    File reference
    File reference_fai
    File QCannotation
  }

  Int space_needed_gb = 10 + round(size([vcf, bam, reference, reference_fai, QCannotation], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    cpu: 1
    memory: "16GB"
    bootDiskSizeGb: 20
    docker: "mgibio/biscuit:0.3.8.2"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /bin/bash /opt/biscuit/scripts/Bisulfite_QC_bisulfiteconversion.sh ~{vcf} ~{bam} ~{reference} ~{QCannotation}
    /bin/bash /opt/biscuit/scripts/Bisulfite_QC_mappingsummary.sh ~{vcf} ~{bam} ~{reference} ~{QCannotation}
    /bin/bash /opt/biscuit/scripts/Bisulfite_QC_CpGretentiondistribution.sh ~{vcf} ~{bam} ~{reference} ~{QCannotation}
    /bin/bash /opt/biscuit/scripts/Bisulfite_QC_Coveragestats.sh ~{vcf} ~{bam} ~{reference} ~{QCannotation}
  >>>

  output {
    # Retention Distribution
    File cpg_retention_dist = "CpGRetentionDist.txt"
    # Conversion
    File base_conversion = "totalBaseConversionRate.txt"
    File read_conversion = "totalReadConversionRate.txt"
    File cph_retention = "CpHRetentionByReadPos.txt"
    File cpg_retention = "CpGRetentionByReadPos.txt"
    # Mapping Summary
    File strand_table = "strand_table.txt"
    File mapping_quality = "mapq_table.txt"
    # Coverage Stats
    File bga_bed = "bga.bed"
    File cov_dist = "covdist_table.txt"
    File bga_bed_dup = "bga_dup.bed"
    File dup_report = "dup_report.txt"
    File cpg_bed = "cpg.bed"
    File cov_dist_cpg = "covdist_cpg_table.txt"
    File cpg_dist = "cpg_dist_table.txt"

    Array[File] qc_files = [ base_conversion, read_conversion,
    cph_retention, cpg_retention, strand_table, mapping_quality,
    cpg_retention_dist, bga_bed, cov_dist, bga_bed_dup, dup_report,
    cpg_bed, cov_dist_cpg, cpg_dist ]
  }
}
