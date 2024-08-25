version 1.0

task mutect {
  input {
    File reference
    File reference_fai
    File reference_dict

    File tumor_bam
    File tumor_bam_bai

    File normal_bam
    File normal_bam_bai

    File interval_list
  }

  parameter_meta {
    tumor_bam: { localization_optional: true }
    tumor_bam_bai: { localization_optional: true }
    normal_bam: { localization_optional: true }
    normal_bam_bai: { localization_optional: true }
  }

  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Float bam_size = size([tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai], "GB")
  Int space_needed_gb = 10 + ceil(reference_size + 2*bam_size + size(interval_list, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "broadinstitute/gatk:4.2.3.0"
    memory: "2GB"
    bootDiskSizeGb: space_needed_gb
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String output_vcf = "mutect.filtered.vcf.gz"
  command <<<
    set -o pipefail

    if [[ `curl metadata.google.internal -i 2> /dev/null | grep 'Metadata-Flavor:'` == "Metadata-Flavor: Google"* ]]; then
       # When the BAMs aren't localized in GCP, samtools needs this token to access them via gs:// URLs.
       export GCS_OAUTH_TOKEN=$(curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token?alt=text | grep access_token | cut -d ' ' -f 2)
    fi

    set -o errexit

    NORMAL=$(samtools view -H ~{normal_bam} | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)
    TUMOR=$(samtools view -H ~{tumor_bam} | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)

    /gatk/gatk Mutect2 --java-options "-Xmx1500m" -O mutect.vcf.gz -R ~{reference} -L ~{interval_list} \
      -I ~{tumor_bam} --read-index ~{tumor_bam_bai} -tumor "$TUMOR" \
      -I ~{normal_bam} --read-index ~{normal_bam_bai} -normal "$NORMAL"

    /gatk/gatk FilterMutectCalls -R ~{reference} -V mutect.vcf.gz -O ~{output_vcf} #Running FilterMutectCalls on the output vcf.
  >>>

  output {
    File vcf = output_vcf
    File vcf_tbi = output_vcf + ".tbi"
  }
}

workflow wf {
  input {
    File reference
    File reference_fai
    File reference_dict
    File tumor_bam
    File tumor_bam_bai
    File normal_bam
    File normal_bam_bai
    File interval_list
  }

  call mutect {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    interval_list=interval_list
  }
}
