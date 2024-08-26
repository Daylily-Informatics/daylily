version 1.0

import "subworkflows/sequence_to_trimmed_fastq.wdl" as sttf
import "tools/bam_to_bigwig.wdl" as btb
import "tools/generate_qc_metrics.wdl" as gqm
import "tools/index_bam.wdl" as ib
import "tools/kallisto.wdl" as k
import "tools/mark_duplicates_and_sort.wdl" as mdas
import "tools/samtools_sort.wdl" as ss
import "tools/star_align_fusion.wdl" as saf
import "tools/star_fusion_detect.wdl" as sfd
import "tools/stringtie.wdl" as s
import "tools/transcript_to_gene.wdl" as ttg
import "tools/xenosplit.wdl" as x
import "types.wdl"

workflow rnaseqStarFusionWithXenosplit {
  input {
    File reference
    File reference_fai
    File reference_dict
    Array[SequenceData] unaligned
    Array[String] outsam_attrrg_line
    File graft_star_genome_dir_zip
    String? graft_outfile_name_prefix
    File host_star_genome_dir_zip
    String? host_outfile_name_prefix
    File star_fusion_genome_dir_zip
    File graft_gtf_file
    File host_gtf_file
    File trimming_adapters
    String trimming_adapter_trim_end
    Int trimming_adapter_min_overlap
    Int trimming_max_uncalled
    Int trimming_min_readlength
    File kallisto_index
    File gene_transcript_lookup_table
    String? strand
    File refFlat
    File ribosomal_intervals
    String sample_name
    Boolean? examine_coding_effect
    String? fusioninspector_mode
    Boolean unzip_fastqs = true
  }

  scatter(sequence in unaligned) {
    call sttf.sequenceToTrimmedFastq {
      input:
      unaligned=sequence,
      adapters=trimming_adapters,
      adapter_trim_end=trimming_adapter_trim_end,
      adapter_min_overlap=trimming_adapter_min_overlap,
      max_uncalled=trimming_max_uncalled,
      min_readlength=trimming_min_readlength,
      unzip_fastqs=unzip_fastqs
    }
  }

  call saf.starAlignFusion as graftStarAlignFusion {
    input:
    outsam_attrrg_line=outsam_attrrg_line,
    star_genome_dir_zip=graft_star_genome_dir_zip,
    outfile_name_prefix=graft_outfile_name_prefix,
    reference_annotation=graft_gtf_file,
    fastq=sequenceToTrimmedFastq.fastq1,
    fastq2=sequenceToTrimmedFastq.fastq2
  }

  call saf.starAlignFusion as hostStarAlignFusion {
    input:
    outsam_attrrg_line=outsam_attrrg_line,
    star_genome_dir_zip=host_star_genome_dir_zip,
    outfile_name_prefix=host_outfile_name_prefix,
    reference_annotation=host_gtf_file,
    fastq=sequenceToTrimmedFastq.fastq1,
    fastq2=sequenceToTrimmedFastq.fastq2
  }

  call x.xenosplit {
    input:
    graftbam=graftStarAlignFusion.aligned_bam,
    hostbam=hostStarAlignFusion.aligned_bam
  }

  call sttf.sequenceToTrimmedFastq as graftbamToFastq {
    input:
    unaligned={"sequence": {"bam": xenosplit.graft_bam}},
    adapters=trimming_adapters,
    adapter_trim_end=trimming_adapter_trim_end,
    adapter_min_overlap=trimming_adapter_min_overlap,
    max_uncalled=trimming_max_uncalled,
    min_readlength=trimming_min_readlength
  }

  call saf.starAlignFusion as graftbamStarAlignFusion {
    input:
    outsam_attrrg_line=outsam_attrrg_line,
    star_genome_dir_zip=graft_star_genome_dir_zip,
    outfile_name_prefix=graft_outfile_name_prefix,
    reference_annotation=graft_gtf_file,
    fastq=[graftbamToFastq.fastq1],
    fastq2=[graftbamToFastq.fastq2]
  }

  call sfd.starFusionDetect {
    input:
    star_fusion_genome_dir_zip=star_fusion_genome_dir_zip,
    examine_coding_effect=examine_coding_effect,
    fusioninspector_mode=fusioninspector_mode,
    fastq=sequenceToTrimmedFastq.fastq1,
    fastq2=sequenceToTrimmedFastq.fastq2,
    outsam_attrrg_line=outsam_attrrg_line
  }

  call k.kallisto {
    input:
    kallisto_index=kallisto_index,
    strand=strand,
    fastqs=[graftbamToFastq.fastqs]
  }

  call ttg.transcriptToGene {
    input:
    transcript_table_h5=kallisto.expression_transcript_h5,
    gene_transcript_lookup_table=gene_transcript_lookup_table
  }

  call ss.samtoolsSort as sortBam {
    input: input_bam=xenosplit.graft_bam
  }

  call mdas.markDuplicatesAndSort as markDup {
    input:
    bam=sortBam.sorted_bam
  }

  call ib.indexBam {
    input: bam=markDup.sorted_bam
  }

  call s.stringtie {
    input:
    bam=indexBam.indexed_bam,
    reference_annotation=graft_gtf_file,
    sample_name=sample_name,
    strand=strand
  }

  call gqm.generateQcMetrics {
    input:
    refFlat=refFlat,
    ribosomal_intervals=ribosomal_intervals,
    strand=strand,
    bam=indexBam.indexed_bam
  }

  call btb.bamToBigwig as cgpbigwigBamcoverage {
    input:
    bam=indexBam.indexed_bam,
    bam_bai=indexBam.indexed_bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict
  }

  output {
    File final_bam = indexBam.indexed_bam
    File final_bam_bai = indexBam.indexed_bam_bai
    File star_fusion_out = graftStarAlignFusion.chim_junc
    File star_junction_out = graftStarAlignFusion.splice_junction_out
    File star_fusion_log = graftStarAlignFusion.log_final
    File star_fusion_predict = starFusionDetect.fusion_predictions
    File star_fusion_abridge = starFusionDetect.fusion_abridged
    File stringtie_transcript_gtf = stringtie.transcript_gtf
    File stringtie_gene_expression_tsv = stringtie.gene_expression_tsv
    File kallisto_transcript_abundance_tsv = kallisto.expression_transcript_table
    File kallisto_transcript_abundance_h5 = kallisto.expression_transcript_h5
    File kallisto_gene_abundance = transcriptToGene.gene_abundance
    File kallisto_fusion_evidence = kallisto.fusion_evidence
    File metrics = generateQcMetrics.metrics
    File? chart = generateQcMetrics.chart
    File xenosplit_statistics = xenosplit.xenosplit_statistics
    File bamcoverage_bigwig = cgpbigwigBamcoverage.outfile
    File? star_fusion_coding_region_effects = starFusionDetect.coding_region_effects
    Array[File] fusioninspector_evidence = starFusionDetect.fusioninspector_evidence
  }
}
