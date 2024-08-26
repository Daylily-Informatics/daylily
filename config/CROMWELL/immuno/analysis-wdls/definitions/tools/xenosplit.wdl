version 1.0

task xenosplit {
  input {
    File graftbam
    File hostbam
  }

  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "20GB"
    bootDiskSizeGb: 100
    docker: "mgibio/xenosplit:0.5"
  }

  command <<<
    set -o pipefail
    set -o errexit

    # Filtering the bam files and preparing them for xenosplit
    /opt/samtools/bin/samtools view -h -F 256 -F 2048 ~{graftbam} | /opt/samtools/bin/samtools sort -n -o graftbam_accepted.bam
    /opt/samtools/bin/samtools view -h -F 256 -F 2048 ~{hostbam} | /opt/samtools/bin/samtools sort -n -o hostbam_accepted.bam

    # Running xenosplit
    python /opt/Xenosplit.py --pairedEnd --out graftOut.bam graftbam_accepted.bam hostbam_accepted.bam
    python /opt/Xenosplit.py --count graftbam_accepted.bam hostbam_accepted.bam > xenosplit_statistics.txt
  >>>

  output {
    File graft_bam = "graftOut.bam"
    File xenosplit_statistics = "xenosplit_statistics.txt"
  }
}
