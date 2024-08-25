version 1.0

task starAlignFusion {
  input {
    Array[File] fastq
    Array[File] fastq2
    # START changed from cwl
    File star_genome_dir_zip
    # END changed from cwl
    Array[String] out_samtype = ["BAM", "Unsorted"]
    String run_mode = "alignReads"
    String out_reads_unmapped = "None"
    Int chim_segment_min = 12
    Int chim_junction_overhang_min = 12
    Int align_sjdb_overhang_min = 10
    Int align_mates_gapmax = 100000
    Int align_intron_max = 100000
    Int chim_segment_read_gapmax = 3
    Array[Int] align_sjstitch_mismatch_nmax = [5, -1, 5, 5]
    String outsam_strand_field = "intronMotif"
    String outsam_unmapped = "Within"
    Array[String] outsam_attrrg_line
    Int chim_multimap_nmax = 10
    Int chim_nonchim_scoredrop_min = 10
    Int peoverlap_nbases_min = 12
    Float peoverlap_mmp = 0.1
    Int chimout_junction_format = 1
    String twopass_mode = "Basic"
    File? reference_annotation
    String outfile_name_prefix = "STAR_"
    String read_files_command = "cat"
    Array[String] outsam_attributes = ["NH", "HI", "AS", "NM", "MD"]
  }

  Int cores = 10
  Float zip_size = size(star_genome_dir_zip, "GB")
  Float fastq_size = size(flatten([fastq, fastq2]), "GB")
  Int space_needed_gb = 10 + round(2 * (zip_size + fastq_size))
  runtime {
    preemptible: 1
    maxRetries: 2
    cpu: cores
    memory: "64GB"
    docker: "trinityctat/starfusion:1.10.1"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String genome_dir = basename(star_genome_dir_zip, ".zip")
  command <<<
    mkdir ~{genome_dir} && unzip -qq ~{star_genome_dir_zip} -d ~{genome_dir}
    /usr/local/bin/STAR --runThreadN ~{cores} \
    --runMode ~{run_mode} \
    --outSAMtype ~{sep=" " out_samtype} \
    --readFilesIn ~{sep="," fastq}  \
    ~{sep="," fastq2}  \
    --outReadsUnmapped ~{out_reads_unmapped} \
    --chimSegmentMin ~{chim_segment_min} \
    --chimJunctionOverhangMin ~{chim_junction_overhang_min} \
    --alignSJDBoverhangMin ~{align_sjdb_overhang_min} \
    --alignMatesGapMax ~{align_mates_gapmax} \
    --alignIntronMax ~{align_intron_max} \
    --chimSegmentReadGapMax ~{chim_segment_read_gapmax} \
    --alignSJstitchMismatchNmax ~{sep=" " align_sjstitch_mismatch_nmax} \
    --outSAMstrandField ~{outsam_strand_field} \
    --outSAMunmapped ~{outsam_unmapped} \
    --outSAMattrRGline ~{sep=" , " outsam_attrrg_line} \
    --chimMultimapNmax ~{chim_multimap_nmax} \
    --chimNonchimScoreDropMin ~{chim_nonchim_scoredrop_min} \
    --peOverlapNbasesMin ~{peoverlap_nbases_min} \
    --peOverlapMMp ~{peoverlap_mmp} \
    --chimOutJunctionFormat ~{chimout_junction_format} \
    --genomeDir ~{genome_dir} \
    --twopassMode ~{twopass_mode} \
    ~{if defined(reference_annotation) then "--sjdbGTFfile ~{select_first([reference_annotation])}" else ""} \
    --outFileNamePrefix ~{outfile_name_prefix} \
    --readFilesCommand ~{read_files_command} \
    --outSAMattributes ~{sep=" " outsam_attributes}
  >>>

  output {
    File aligned_bam = "~{outfile_name_prefix}Aligned.out.bam"
    File log_final = "~{outfile_name_prefix}Log.final.out"
    File log = "~{outfile_name_prefix}Log.out"
    File log_progress = "~{outfile_name_prefix}Log.progress.out"
    File splice_junction_out = "~{outfile_name_prefix}SJ.out.tab"
    File chim_junc = "~{outfile_name_prefix}Chimeric.out.junction"
  }
}

workflow wf { call starAlignFusion { input: } }
