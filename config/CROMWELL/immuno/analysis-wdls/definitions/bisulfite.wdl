version 1.0

import "subworkflows/bam_to_trimmed_fastq_and_biscuit_alignments.wdl" as bttfaba
import "tools/bisulfite_qc.wdl" as bq
import "tools/bam_to_cram.wdl" as btc
import "tools/bedgraph_to_bigwig.wdl" as btb
import "tools/biscuit_pileup.wdl" as bp
import "tools/bisulfite_vcf2bed.wdl" as bv
import "tools/index_cram.wdl" as ic
import "tools/merge_bams.wdl" as mb

workflow bisulfite {
  input {
    File reference
    File reference_fai
    File reference_dict
    File reference_sizes
    Array[File] instrument_data_bams
    Array[String] read_group_id
    File trimming_adapters
    String trimming_adapter_trim_end
    Int trimming_adapter_min_overlap
    Int trimming_max_uncalled
    Int trimming_min_readlength
    File QCannotation
    Boolean assay_non_cpg_sites = false
  }

  scatter(bam in instrument_data_bams) {
    scatter(rgi in read_group_id) {
      call bttfaba.bamToTrimmedFastqAndBiscuitAlignments {
        input:
        bam=bam,
        read_group_id=rgi,
        adapters=trimming_adapters,
        adapter_trim_end=trimming_adapter_trim_end,
        adapter_min_overlap=trimming_adapter_min_overlap,
        max_uncalled=trimming_max_uncalled,
        min_readlength=trimming_min_readlength,
        reference_index=reference
      }
    }
  }

  call mb.mergeBams as merge {
    input:
    bams=flatten(bamToTrimmedFastqAndBiscuitAlignments.aligned_bam)
  }

  call bp.biscuitPileup as pileup {
    input:
    bam=merge.merged_bam,
    reference=reference,
    reference_fai=reference_fai
  }

  call bq.bisulfiteQc {
    input:
    vcf=pileup.vcf,
    bam=merge.merged_bam,
    reference=reference,
    reference_fai=reference_fai,
    QCannotation=QCannotation
  }

  call bv.bisulfiteVcf2bed as vcf2bed {
    input:
    vcf=pileup.vcf,
    reference=reference,
    reference_fai=reference_fai,
    assay_non_cpg_sites=assay_non_cpg_sites
  }

  call btb.bedgraphToBigwig {
    input:
    methylation_bedgraph=vcf2bed.methylation_bedgraph,
    reference_sizes=reference_sizes
  }

  call btc.bamToCram {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=merge.merged_bam
  }

  call ic.indexCram {
    input:
    cram=bamToCram.cram
  }

  output {
    File cram = indexCram.indexed_cram
    File cram_crai = indexCram.indexed_cram_crai
    File vcf = pileup.vcf
    Array[File] cpgs = vcf2bed.methylation_bed
    Array[File] cpg_bigwig = bedgraphToBigwig.methylation_bigwig
    Array[File] qc_files = bisulfiteQc.qc_files
  }
}
