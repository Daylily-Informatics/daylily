version 1.0

task samtoolsFlagstat {
  input {
    File bam
    File bam_bai
  }

  parameter_meta {
    bam: { localization_optional: true }
    bam_bai: { localization_optional: true }
  }

  Int space_needed_gb = 10 + round(size([bam, bam_bai], "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = basename(bam) + ".flagstat"
  command <<<
    ACCESS_TOKEN=$(wget -O - --header "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token?alt=text 2> /dev/null | grep access_token)
    if [[ "$ACCESS_TOKEN" == "access_token"* ]]; then
       # When the BAMs aren't localized in GCP, samtools needs this token to access them via gs:// URLs.
       export GCS_OAUTH_TOKEN=$(echo "$ACCESS_TOKEN" | cut -d ' ' -f 2 )
       echo "got token" ${GCS_OAUTH_TOKEN:0:5}
    fi

    /usr/local/bin/samtools flagstat ~{bam} > ~{outfile}
  >>>

  output {
    File flagstats = outfile
  }
}

workflow wf {
  input {
    File bam
    File bam_bai
  }

  call samtoolsFlagstat {
    input:
    bam=bam,
    bam_bai=bam_bai
  }
}
