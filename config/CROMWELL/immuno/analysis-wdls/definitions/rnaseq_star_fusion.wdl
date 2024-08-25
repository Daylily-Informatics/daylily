version 1.0

import "subworkflows/sequence_to_trimmed_fastq.wdl" as sttf
import "tools/agfusion.wdl" as a
import "tools/bam_to_cram.wdl" as btc
import "tools/generate_qc_metrics.wdl" as gqm
import "tools/index_bam.wdl" as ib
import "tools/index_cram.wdl" as ic
import "tools/kallisto.wdl" as k
import "tools/mark_duplicates_and_sort.wdl" as mdas
import "tools/samtools_sort.wdl" as ss
import "tools/star_fusion_detect.wdl" as sfd
import "tools/strandedness_check.wdl" as sc
import "tools/stringtie.wdl" as s
import "tools/transcript_to_gene.wdl" as ttg
import "tools/samtools_flagstat.wdl" as sf
import "types.wdl"  # !UnusedImport

workflow rnaseqStarFusion {
  input {
    File reference
    File reference_fai
    File reference_dict

    Array[SequenceData] unaligned
    String? strand  # enum [first, second, unstranded]
    String sample_name
    
    File star_fusion_genome_dir_zip
    File cdna_fasta
    File reference_annotation

    File trimming_adapters
    String trimming_adapter_trim_end
    Int trimming_adapter_min_overlap
    Int trimming_max_uncalled
    Int trimming_min_readlength

    File kallisto_index
    File gene_transcript_lookup_table
    File refFlat
    File ribosomal_intervals
    Boolean unzip_fastqs = true

    Boolean? examine_coding_effect
    String? fusioninspector_mode  # enum [inspect validate]
    File agfusion_database
    Boolean? agfusion_annotate_noncanonical
    Float? min_ffpm_level
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

    call sc.strandednessCheck {
      input:
      reference_annotation=reference_annotation,
      kallisto_index=kallisto_index,
      cdna_fasta=cdna_fasta,
      reads1=sequenceToTrimmedFastq.fastq1,
      reads2=sequenceToTrimmedFastq.fastq2
    }

    String? attrrg_line = sequence.readgroup
  }

  call sfd.starFusionDetect {
    input:
    star_fusion_genome_dir_zip=star_fusion_genome_dir_zip,
    examine_coding_effect=examine_coding_effect,
    fusioninspector_mode=fusioninspector_mode,
    fastq=sequenceToTrimmedFastq.fastq1,
    fastq2=sequenceToTrimmedFastq.fastq2,
    outsam_attrrg_line=select_all(attrrg_line),
    min_ffpm_level=min_ffpm_level
  }

  call k.kallisto {
    input:
    kallisto_index=kallisto_index,
    strand=strand,
    fastqs=sequenceToTrimmedFastq.fastqs
  }

  call ttg.transcriptToGene {
    input:
    transcript_table_h5=kallisto.expression_transcript_h5,
    gene_transcript_lookup_table=gene_transcript_lookup_table
  }

  call ss.samtoolsSort as sortBam {
    input: input_bam=starFusionDetect.aligned_bam
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
    bam=markDup.sorted_bam,
    reference_annotation=reference_annotation,
    sample_name=sample_name,
    strand=strand
  }

  call gqm.generateQcMetrics {
    input:
    refFlat=refFlat,
    ribosomal_intervals=ribosomal_intervals,
    strand=strand,
    bam=markDup.sorted_bam
  }

  call btc.bamToCram {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=indexBam.indexed_bam
  }

  call ic.indexCram {
    input: cram=bamToCram.cram
  }

  call a.agfusion {
    input:
    fusion_predictions=starFusionDetect.fusion_predictions,
    agfusion_database=agfusion_database,
    annotate_noncanonical=agfusion_annotate_noncanonical
  }

  call sf.samtoolsFlagstat {
    input:
    bam=indexBam.indexed_bam,
    bam_bai=indexBam.indexed_bam_bai
  }

  output {
    File cram = indexCram.indexed_cram
    File cram_crai = indexCram.indexed_cram_crai
    File star_fusion_out = starFusionDetect.chim_junc
    File star_junction_out = starFusionDetect.splice_junction_out
    File star_fusion_log_final = starFusionDetect.log_final
    File star_fusion_log = starFusionDetect.log
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
    Array[File] strand_info = strandednessCheck.strandedness_check
    File final_bam = indexBam.indexed_bam
    File final_bam_bai = indexBam.indexed_bam_bai
    File final_bai = indexBam.indexed_bai
    File annotated_fusion_predictions_zip = agfusion.annotated_fusion_predictions_zip
    File? star_fusion_coding_region_effects = starFusionDetect.coding_region_effects
    Array[File] fusioninspector_evidence = starFusionDetect.fusioninspector_evidence
    File flagstats = samtoolsFlagstat.flagstats
    PreliminaryStarFusionResults prelim_starfusion_results = starFusionDetect.prelim_results
  }
}
