version 1.0

task bamToCram {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
  }

  parameter_meta {
    bam: { localization_optional: true }
  }

  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Int size_needed_gb = 10 + round(size(bam, "GB") * 2 + reference_size)
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    memory: "4GB"
    disks: "local-disk ~{size_needed_gb} HDD"
  }

  String outfile = basename(bam, ".bam") + ".cram"
  command <<<
    ACCESS_TOKEN=$(wget -O - --header "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token?alt=text 2> /dev/null | grep access_token)
    if [[ "$ACCESS_TOKEN" == "access_token"* ]]; then
       # When the BAMs aren't localized in GCP, samtools needs this token to access them via gs:// URLs.
       export GCS_OAUTH_TOKEN=$(echo "$ACCESS_TOKEN" | cut -d ' ' -f 2 )
       echo "got token" ${GCS_OAUTH_TOKEN:0:5}
    fi

    /usr/local/bin/samtools view -C -T ~{reference} ~{bam} > ~{outfile}
  >>>

  output {
    File cram = outfile
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict
    File bam
  }

  call bamToCram {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    bam=bam
  }
}
