version 1.0

import "tools/downsample.wdl" as d
import "subworkflows/gatk_haplotypecaller_iterator.wdl" as ghi
import "tools/collect_wgs_metrics.wdl" as cwm

struct Downsample {
  File cram
  Float downsample_ratio
  Float contamination
}

workflow downsampleAndRecall {
  input {
    File reference
    File reference_fai
    File reference_dict
    Array[Downsample] crams_to_downsample
    String? downsample_strategy  # enum [HighAccuracy, ConstantMemory, Chained]
    Int? downsample_seed
    String emit_reference_confidence  # enum [NONE, BP_RESOLUTION, GVCF]
    Int? max_alternate_alleles
    Int? ploidy
    String? read_filter
    Array[Array[String]] intervals
    Int qc_minimum_mapping_quality
    Int qc_minimum_base_quality
  }

  scatter(x in crams_to_downsample) {
    call d.downsample {
      input:
      sam=x.cram,
      probability=x.downsample_ratio,
      reference=reference,
      reference_fai=reference_fai,
      reference_dict=reference_dict,
      random_seed=downsample_seed,
      strategy=downsample_strategy
    }

    call ghi.gatkHaplotypecallerIterator as haplotypeCaller {
      input:
      reference=reference,
      reference_fai=reference_fai,
      reference_dict=reference_dict,
      bam=downsample.downsampled_sam,
      bai=downsample.downsampled_sam_bai,
      emit_reference_confidence=emit_reference_confidence,
      gvcf_gq_bands=[],
      intervals=intervals,
      contamination_fraction=write_lines(["~{x.contamination}"]),
      max_alternate_alleles=max_alternate_alleles,
      ploidy=ploidy,
      read_filter=read_filter,
      output_prefix=basename(downsample.downsampled_sam, ".bam") + ".downsampled."
    }

    call cwm.collectWgsMetrics {
      input:
      bam=downsample.downsampled_sam,
      bam_bai=downsample.downsampled_sam_bai,
      reference=reference,
      reference_fai=reference_fai,
      reference_dict=reference_dict,
      minimum_mapping_quality=qc_minimum_mapping_quality,
      minimum_base_quality=qc_minimum_base_quality,
      sample_name=basename(downsample.downsampled_sam, ".bam")
    }
  }

  output {
    Array[Array[File]] gvcfs = haplotypeCaller.gvcf
    Array[File] wgs_metrics = collectWgsMetrics.wgs_metrics
  }
}
