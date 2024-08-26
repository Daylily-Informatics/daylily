version 1.0

struct PreliminaryStarFusionResults {
  File candidates
  File bp_filter
  File candidates_filtered_ffpm
  File wannot
  File wannot_filter
  File wannot_filter_artifact
  File wannot_filter_artifact_minffpm
}

task starFusionDetect {
  input {
    File star_fusion_genome_dir_zip
    String fusion_output_dir = "STAR-Fusion_outdir"
    String star_path = "/usr/local/bin/STAR"
    # TODO: is this presence or =true ?
    Boolean examine_coding_effect = false
    String? fusioninspector_mode  # enum [inspect, validate]
    Array[File] fastq
    Array[File] fastq2
    Array[String] outsam_attrrg_line
    Float? min_ffpm_level = 0.05
  }

  Int cores = 24
  Float zip_size = size(star_fusion_genome_dir_zip, "GB")
  Float fastq_size = size(flatten([fastq, fastq2]), "GB")
  Int space_needed_gb = 10 + round(3 * (zip_size + fastq_size))
  runtime {
    preemptible: 0
    maxRetries: 2
    memory: "64GB"
    cpu: cores
    docker: "mgibio/starfusion:pre-1.11.c-fefaf71"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

    # https://github.com/STAR-Fusion/STAR-Fusion/issues/175#issuecomment-567913451
  String genome_lib_dir = "`pwd`/" + basename(star_fusion_genome_dir_zip, ".zip")
  command <<<
    mkdir ~{genome_lib_dir} && unzip -qq ~{star_fusion_genome_dir_zip} -d ~{genome_lib_dir}
    /usr/local/src/STAR-Fusion/STAR-Fusion --CPU ~{cores} \
        --genome_lib_dir ~{genome_lib_dir} \
        --output_dir ~{fusion_output_dir} --STAR_PATH ~{star_path} \
        ~{true="--examine_coding_effect" false="" examine_coding_effect} \
        ~{if defined(fusioninspector_mode) then "--FusionInspector " + fusioninspector_mode else ""} \
        --STAR_outSAMattrRGline "~{sep=" , " outsam_attrrg_line}" \
        --left_fq ~{sep="," fastq} --right_fq ~{sep="," fastq2} \
        --min_FFPM ~{min_ffpm_level}
  >>>

  output {
    # star fusion outputs 
    File fusion_predictions = fusion_output_dir + "/star-fusion.fusion_predictions.tsv"
    File fusion_abridged = fusion_output_dir + "/star-fusion.fusion_predictions.abridged.tsv"
    File? coding_region_effects = fusion_output_dir + "/star-fusion.fusion_predictions.abridged.coding_effect.tsv"
    # Fusion inspector outputs
    # if no mode specified, this will just not find any files
    Array[File] fusioninspector_evidence = glob(fusion_output_dir + "/FusionInspector-" + select_first([fusioninspector_mode, ""]) + "/finspector.*")
    File fusioninspector_log = fusion_output_dir + "/FusionInspector.log"
    # STAR alignment files
    File aligned_bam = fusion_output_dir + "/Aligned.out.bam"
    File log_final = fusion_output_dir + "/Log.final.out"
    File log = fusion_output_dir + "/Log.out"
    File log_progress = fusion_output_dir + "/Log.progress.out"
    File splice_junction_out = fusion_output_dir + "/SJ.out.tab"
    File chim_junc = fusion_output_dir + "/Chimeric.out.junction"
    # STAR also outputs gene counts file just like Kallisto
    File gene_counts = fusion_output_dir + "/ReadsPerGene.out.tab"

    PreliminaryStarFusionResults prelim_results = object {
      candidates: fusion_output_dir + "/star-fusion.preliminary/star-fusion.fusion_candidates.preliminary",
      bp_filter: fusion_output_dir + "/star-fusion.preliminary/star-fusion.filter.intermediates_dir/star-fusion.post_blast_and_promiscuity_filter",
      candidates_filtered_ffpm: fusion_output_dir + "/star-fusion.preliminary/star-fusion.fusion_candidates.preliminary.filtered.FFPM",
      wannot: fusion_output_dir + "/star-fusion.preliminary/star-fusion.fusion_candidates.preliminary.wSpliceInfo.wAnnot",
      wannot_filter: fusion_output_dir + "/star-fusion.preliminary/star-fusion.fusion_candidates.preliminary.wSpliceInfo.wAnnot.annot_filter.pass",
      wannot_filter_artifact: fusion_output_dir + "/star-fusion.preliminary/star-fusion.fusion_candidates.preliminary.wSpliceInfo.wAnnot.annot_filter.pass.RTartifact.pass",
      wannot_filter_artifact_minffpm: fusion_output_dir + "/star-fusion.preliminary/star-fusion.fusion_candidates.preliminary.wSpliceInfo.wAnnot.annot_filter.pass.RTartifact.pass.minFFPM.${min_ffpm_level}.pass"
    }
  }
}

workflow wf { call starFusionDetect { input: } }
