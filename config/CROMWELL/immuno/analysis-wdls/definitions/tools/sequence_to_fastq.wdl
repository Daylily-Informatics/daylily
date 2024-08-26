version 1.0

import "../types.wdl"  # !UnusedImport

task sequenceToFastq {
  input {
    File? bam
    File? fastq1
    File? fastq2
    Boolean unzip_fastqs = false
  }

  Int compression_multiplier = if unzip_fastqs then 10 else 1
  Int space_needed_gb = 10 + ceil(2*compression_multiplier*size([bam, fastq1, fastq2], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "16GB"
    bootDiskSizeGb: 25
    docker: "broadinstitute/picard:2.23.6"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outdir = "outdir"
  command <<<
    set -o pipefail
    set -o errexit
    set -o nounset

    UNZIP=~{unzip_fastqs}
    OUTDIR=~{outdir}
    BAM=~{if defined(bam) then bam else ""}
    FASTQ1=~{if defined(fastq1) then fastq1 else ""}
    FASTQ2=~{if defined(fastq2) then fastq2 else ""}
    MODE=~{if defined(bam) then "bam" else "fastq"}

    mkdir -p $OUTDIR

    if [[ "$MODE" == "fastq" ]]; then #must be fastq input
        if $UNZIP; then
            if gzip -t $FASTQ1 2> /dev/null; then
                gunzip -c $FASTQ1 > $OUTDIR/read1.fastq
            else
                cp $FASTQ1 $OUTDIR/read1.fastq
            fi

            if gzip -t $FASTQ2 2> /dev/null; then
                gunzip -c $FASTQ2 > $OUTDIR/read2.fastq
            else
                cp $FASTQ2 $OUTDIR/read2.fastq
            fi
        else
            if gzip -t $FASTQ1 2> /dev/null; then
                cp $FASTQ1 $OUTDIR/read1.fastq.gz
            else
                cp $FASTQ1 $OUTDIR/read1.fastq
            fi

            if gzip -t $FASTQ2 2> /dev/null; then
                cp $FASTQ2 $OUTDIR/read2.fastq.gz
            else
                cp $FASTQ2 $OUTDIR/read2.fastq
            fi
        fi

    else # then
        ##run samtofastq here, dumping to the same filenames
        ## input file is $BAM
        /usr/bin/java -Xmx4g -jar /usr/picard/picard.jar SamToFastq I="$BAM" INCLUDE_NON_PF_READS=true F=$OUTDIR/read1.fastq F2=$OUTDIR/read2.fastq VALIDATION_STRINGENCY=SILENT
    fi
  >>>

  output {
    File read1_fastq = glob("~{outdir}/read1.fastq*")[0]
    File read2_fastq = glob("~{outdir}/read2.fastq*")[0]
  }
}

workflow wf {
  input {
    File? bam
    File? fastq1
    File? fastq2
    Boolean unzip_fastqs = false
  }

  call sequenceToFastq {
    input:
    bam=bam,
    fastq1=fastq1,
    fastq2=fastq2,
    unzip_fastqs=unzip_fastqs
  }
}
