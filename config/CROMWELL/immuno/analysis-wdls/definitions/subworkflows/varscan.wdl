version 1.0

import "../tools/varscan_somatic.wdl" as vs
import "../tools/varscan_process_somatic.wdl" as vps

workflow varscan {
  input {
    File reference
    File reference_fai
    File reference_dict

    File tumor_bam
    File tumor_bam_bai

    File normal_bam
    File normal_bam_bai

    File? roi_bed

    Int? strand_filter
    Int? min_coverage
    Float? varscan_min_var_freq
    Float? p_value
    Float? max_normal_freq
  }

  call vs.varscanSomatic as somatic {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    roi_bed=roi_bed,
    strand_filter=strand_filter,
    min_coverage=min_coverage,
    varscan_min_var_freq=varscan_min_var_freq,
    p_value=p_value
  }

  call vps.varscanProcessSomatic as process_somatic_snvs {
    input:
    variants=somatic.snvs,
    max_normal_freq=max_normal_freq
  }

  call vps.varscanProcessSomatic as process_somatic_indels {
    input:
    variants=somatic.indels,
    max_normal_freq=max_normal_freq
  }

  output {
    # somatic
    File snvs = somatic.snvs
    File indels = somatic.indels
    # process_somatic_snvs
    File somatic_hc_snvs = process_somatic_snvs.somatic_hc
    File somatic_snvs = process_somatic_snvs.somatic
    File germline_hc_snvs = process_somatic_snvs.germline_hc
    File germline_snvs = process_somatic_snvs.germline
    File loh_hc_snvs = process_somatic_snvs.loh_hc
    File loh_snvs = process_somatic_snvs.loh
    # process_somatic_indels
    File somatic_hc_indels = process_somatic_indels.somatic_hc
    File somatic_indels = process_somatic_indels.somatic
    File germline_hc_indels = process_somatic_indels.germline_hc
    File germline_indels = process_somatic_indels.germline
    File loh_hc_indels = process_somatic_indels.loh_hc
    File loh_indels = process_somatic_indels.loh
  }
}
