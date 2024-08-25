version 1.0

import "../types.wdl"  # !UnusedImport

task sequenceAlignAndTag {
  input {
    SequenceData unaligned
    TrimmingOptions? trimming
    File reference
    File reference_alt
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_0123
  }

  # Disk space
  Float data_size = size([unaligned.sequence.bam, unaligned.sequence.fastq1, unaligned.sequence.fastq2], "GB")
  Float reference_size = size([reference, reference_alt, reference_amb, reference_ann, reference_bwt, reference_pac, reference_0123], "GB")
  Int space_needed_gb = 10 + ceil(5*data_size + reference_size)
  # CPU |  Memory / RAM
  #24 bwa cores + 1 for samblaster + 2 for samtools view
  Int cores = 64 
  # mem values tested on a 60x WGS bam with up to 30 cores. Lower values would 
  # be fine for exome if we ever want to optimize to that extent
  # Mem does not seem to be as dependent on # of cores as I expected
  Int instance_memory_gb = 76 
  Int jvm_memory_gb = 4
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/alignment_helper-cwl:2.2.1"
    memory: "~{instance_memory_gb}GB"
    cpu: cores
    bootDiskSizeGb: space_needed_gb
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outname = "refAligned.bam"
  command <<<
    set -o pipefail
    set -o errexit
    set -o nounset

    # destructure unaligned
    MODE=~{if defined(unaligned.sequence.bam) then "bam" else "fastq" }
    BAM="~{unaligned.sequence.bam}"
    FASTQ1="~{unaligned.sequence.fastq1}"
    FASTQ2="~{unaligned.sequence.fastq2}"
    REFERENCE="~{reference}"
    READGROUP="~{unaligned.readgroup}"
    # destructure trimming
    RUN_TRIMMING=~{if defined(trimming) then "true" else "false"}
    TRIMMING_ADAPTERS=~{if defined(trimming) then "~{select_first([trimming]).adapters}" else ""}
    TRIMMING_ADAPTER_MIN_OVERLAP=~{if defined(trimming) then "~{select_first([trimming]).min_overlap}" else ""}

    function bwa_blast_view () {
        /usr/local/bin/bwa-mem2 mem -K 100000000 -t ~{cores} -Y -p -R "$READGROUP" "$REFERENCE" /dev/stdin \
            | /usr/local/bin/samblaster -a --addMateTags \
            | /opt/samtools/bin/samtools view -@ 2 -b -S /dev/stdin
    }

    if [[ "$MODE" == "fastq" ]]; then
        if [[ "$RUN_TRIMMING" == "false" ]]; then
            /usr/local/bin/bwa-mem2 mem -K 100000000 -t ~{cores} -Y -R "$READGROUP" "$REFERENCE" $FASTQ1 $FASTQ2 \
                | /usr/local/bin/samblaster -a --addMateTags \
                | /opt/samtools/bin/samtools view -@ 2 -b -S /dev/stdin > "~{outname}"
        else
            /opt/flexbar/flexbar --adapters "$TRIMMING_ADAPTERS" --reads $FASTQ2 --reads2 $FASTQ2 --adapter-trim-end LTAIL --adapter-min-overlap "$TRIMMING_ADAPTER_MIN_OVERLAP" --adapter-error-rate 0.1 --max-uncalled 300 --stdout-reads \
                | bwa_blast_view > "~{outname}"
        fi
    fi
    if [[ "$MODE" == "bam" ]]; then
        if [[ "$RUN_TRIMMING" == "false" ]]; then
            /usr/bin/java -Xmx~{jvm_memory_gb}g -jar /opt/picard/picard.jar SamToFastq I=$BAM INTERLEAVE=true INCLUDE_NON_PF_READS=true FASTQ=/dev/stdout \
                | bwa_blast_view > "~{outname}"
        else
            /usr/bin/java -Xmx~{jvm_memory_gb}g -jar /opt/picard/picard.jar SamToFastq I=$BAM INTERLEAVE=true INCLUDE_NON_PF_READS=true FASTQ=/dev/stdout \
                | /opt/flexbar/flexbar --adapters "$TRIMMING_ADAPTERS" --reads - --interleaved --adapter-trim-end LTAIL --adapter-min-overlap "$TRIMMING_ADAPTER_MIN_OVERLAP" --adapter-error-rate 0.1 --max-uncalled 300 --stdout-reads \
                | bwa_blast_view > "~{outname}"
        fi
    fi
  >>>

  output {
    File aligned_bam = outname
  }
}

workflow wf {
  input {
    SequenceData unaligned
    TrimmingOptions? trimming
    File reference
    File reference_alt
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_0123
  }
  call sequenceAlignAndTag {
    input:
    unaligned=unaligned,
    trimming=trimming,
    reference=reference,
    reference_alt=reference_alt,
    reference_amb=reference_amb,
    reference_ann=reference_ann,
    reference_bwt=reference_bwt,
    reference_pac=reference_pac,
    reference_0123=reference_0123
  }
}
