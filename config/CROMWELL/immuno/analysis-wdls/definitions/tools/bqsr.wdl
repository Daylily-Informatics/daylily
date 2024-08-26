version 1.0

# takes an aligned, sorted bam and returns a bam to which base quality
# score recalibration has been applied

workflow doBqsr {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bam_bai
    Array[File] known_sites
    Array[File] known_sites_tbi
    String output_name = "final"
    Int preemptible_tries = 3
  }

  Float bam_size = size([bam, bam_bai], "GB")
  # Create list of sequences for scattering
  call CreateSequenceGroupingTSV {
    input:
      reference_dict = reference_dict,
      preemptible_tries = preemptible_tries
  }

  #create the bqsr model, scattering over intervals
  scatter (interval in CreateSequenceGroupingTSV.sequence_grouping) {
    call bqsr {
      input:
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        bam = bam,
        bam_bai = bam_bai,
        known_sites = known_sites,
        known_sites_tbi = known_sites_tbi,
        intervals = interval,
        preemptible_tries = preemptible_tries
    }
  }

  # Merge the bqsr model reports from the scatter
  call GatherBqsrReports {
    input:
      input_bqsr_reports = bqsr.bqsr_table,
      preemptible_tries = preemptible_tries
  }

  scatter (interval in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
    call applyBqsr {
      input:
        reference=reference,
        reference_fai=reference_fai,
        reference_dict=reference_dict,
        bam=bam,
        bam_bai=bam_bai,
        bqsr_table=GatherBqsrReports.bqsr_table,
        output_name=output_name,
        intervals = interval,
        preemptible_tries = preemptible_tries
    }
  }

  call GatherBamFiles {
    input:
      input_bams = applyBqsr.bqsr_bam,
      output_bam_basename = output_name,
      bam_size = bam_size,
      preemptible_tries = preemptible_tries,
  }

  output {
    File bqsr_bam = GatherBamFiles.output_bam
    File bqsr_bai = GatherBamFiles.output_bai
    File bqsr_bam_bai = GatherBamFiles.output_bam_bai
  }
}

#------------------------------------------------------------

# Generate sets of intervals for scatter-gathering over chromosomes
# Code modified from https://github.com/gatk-workflows/gatk4-data-processing/blob/master/processing-for-variant-discovery-gatk4.wdl
# keeps whole chromosomes intact, just tries to balance the load by grouping them.
# e.g. for GRCh38, creates 18 groups
task CreateSequenceGroupingTSV {
 input {
    File reference_dict
    Int preemptible_tries
  }

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("~{reference_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in
    # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
    #hg38_protection_tag = ":1+"
    hg38_protection_tag = ""
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  runtime {
    maxRetries: 2
    preemptible: preemptible_tries
    docker: "python:2.7"
    memory: "2 GiB"
  }
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}

#run the base recalibrator to generate the metrics table
task bqsr {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bam_bai
    Array[File] known_sites
    Array[File] known_sites_tbi
    Array[String] intervals
    Int preemptible_tries
  }
  Float bam_size = size([bam, bam_bai], "GB")
  Float known_sites_size = size(known_sites, "GB") + size(known_sites_tbi, "GB")
  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Int space_needed_gb = 10 + round(known_sites_size  + bam_size + reference_size)
  runtime {
    maxRetries: 2
    docker: "broadinstitute/gatk:4.1.8.1"
    memory: "6GB"
    disks: "local-disk ~{space_needed_gb} SSD"
    preemptible: preemptible_tries
  }

  String outfile = "bqsr.table"
 command { 
    /gatk/gatk --java-options -Xmx4g BaseRecalibrator \
    -O ~{outfile} \
    -L ~{sep=" -L " intervals} \
    -R ~{reference} \
    -I ~{bam} \
    --known-sites ~{sep=" --known-sites " known_sites} 
  }

  output {
    File bqsr_table = outfile
    Float input_bam_size = bam_size
  }
}

#pull bqsr reports from the scattered runs into a single report
task GatherBqsrReports {
 input {
   Array[File] input_bqsr_reports
   Int preemptible_tries
  }
  command {
    /gatk/gatk --java-options -Xmx4g GatherBQSRReports \
      -I ~{sep=' -I ' input_bqsr_reports} \
      -O bqsr_report.txt
  }
  runtime {
    maxRetries: 2
    preemptible: preemptible_tries
    docker: "broadinstitute/gatk:4.1.8.1"
    memory: "4 GiB"
    disks: "local-disk 10 HDD"
  }
  output {
    File bqsr_table = "bqsr_report.txt"
  }
}

#use the metrics table to update the bam
task applyBqsr {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
    File bam_bai
    File bqsr_table
    String output_name = "bqsr"
    Array[String] intervals
    Int preemptible_tries
  }

  Int space_needed_gb = 10 + round(size([bqsr_table, reference, reference_fai, reference_dict], "GB") + size([bam, bam_bai], "GB") * 2)
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "broadinstitute/gatk:4.1.8.1"
    memory: "18GB"
    disks: "local-disk ~{space_needed_gb} SSD"
  }

  command <<<
    /gatk/gatk --java-options -Xmx16g ApplyBQSR \
    -O ~{output_name}.bam \
    -L ~{sep=" -L " intervals} \
    ~{sep=" " prefix("--static-quantized-quals ", [10, 20, 30])} \
    -R ~{reference} \
    -I ~{bam} \
    -bqsr ~{bqsr_table}
  >>>

  output {
    File bqsr_bam = "~{output_name}.bam"
  }
}


# Combine the scattered bams from applyBqsr
task GatherBamFiles {
  input {
    Array[File] input_bams
    String output_bam_basename = "final"
    Float mem_size_gb = 3
    Float bam_size
    Int preemptible_tries
  }
  Int command_mem_gb = ceil(mem_size_gb) - 1
  Int space_needed_gb = 10 + round(bam_size * 2.5)
  runtime {
    maxRetries: 2
    preemptible: preemptible_tries
    docker: "broadinstitute/gatk:4.1.8.1"
    memory: "3 GiB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }
  command {
     /gatk/gatk --java-options "-Xmx2G" GatherBamFiles \
      --INPUT ~{sep=' --INPUT ' input_bams} \
      --OUTPUT ~{output_bam_basename}.bam \
      --CREATE_INDEX true \
    # if we want md5s in the future, this is the place
    #  --CREATE_MD5_FILE true
    # some tools require file.bai, some file.bam.bai - generate both
    cp ~{output_bam_basename}.bai ~{output_bam_basename}.bam.bai
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_bai = "~{output_bam_basename}.bam.bai"
    File output_bai = "~{output_bam_basename}.bai"
    # File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}
