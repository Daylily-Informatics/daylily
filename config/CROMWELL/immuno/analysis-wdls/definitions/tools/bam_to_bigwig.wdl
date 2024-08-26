version 1.0

task bamToBigwig {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
  }

  Float bam_size_gb = size([bam, bam_bai], "GB")
  Float reference_size_gb = size([reference, reference_fai, reference_dict], "GB")
  Int space_needed_gb = 10 + round(3*bam_size_gb + reference_size_gb)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "32GB"
    docker: "quay.io/biocontainers/cgpbigwig:1.4.0--h93d22ca_0"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String output_bw = basename(bam, ".bam") + ".bw"
  command <<<
    bam2bw -a -F 1024 -o ~{output_bw} -i ~{bam} -r ~{reference}
  >>>

  output {
    File outfile = output_bw
  }
}


workflow wf {
  input {
    File bam
    File bam_bai
    File reference
    File reference_fai
    File reference_dict
  }

  call bamToBigwig {
    input:
    bam=bam,
    bam_bai=bam_bai,
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict
  }
}
