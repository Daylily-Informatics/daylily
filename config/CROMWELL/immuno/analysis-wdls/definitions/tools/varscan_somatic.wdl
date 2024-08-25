version 1.0

task varscanSomatic {
  input {
    File reference
    File reference_fai
    File reference_dict

    File tumor_bam
    File tumor_bam_bai

    File normal_bam
    File normal_bam_bai

    Int strand_filter = 0
    Int min_coverage = 8
    Float varscan_min_var_freq = 0.05
    Float p_value = 0.99
    File? roi_bed
  }

  Float reference_size = size([reference, reference_fai, reference_dict], "GB")
  Float bam_size = size([tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai], "GB")
  Int space_needed_gb = 10 + ceil(reference_size + bam_size*2)
  runtime {
    preemptible: 1
    maxRetries: 2

    memory: "12GB"
    cpu: 2
    docker: "mgibio/varscan-cwl:v2.4.2-samtools1.16.1"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
 
    ACCESS_TOKEN=$(wget -O - --header "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token?alt=text 2> /dev/null | grep access_token)
    if [[ "$ACCESS_TOKEN" == "access_token"* ]]; then
       # When the BAMs aren't localized in GCP, samtools needs this token to access them via gs:// URLs.
       export GCS_OAUTH_TOKEN=$(echo "$ACCESS_TOKEN" | cut -d ' ' -f 2 )
       echo "got token" ${GCS_OAUTH_TOKEN:0:5}
    fi

    set -o errexit
    set -o nounset
    set -o pipefail


    /opt/samtools/bin/samtools mpileup --no-baq ~{if defined(roi_bed) then "-l ~{roi_bed}" else ""} -f "~{reference}" "~{normal_bam}" "~{tumor_bam}" | \
    java -jar /opt/varscan/VarScan.jar somatic /dev/stdin \
    "output" \
    --strand-filter "~{strand_filter}" \
    --min-coverage "~{min_coverage}" \
    --min-var-freq "~{varscan_min_var_freq}" \
    --p-value "~{p_value}" \
    --mpileup 1 \
    --output-vcf
  >>>

  output {
    File snvs = "output.snp.vcf"
    File indels = "output.indel.vcf"
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

    Int strand_filter
    Int min_coverage
    Float varscan_min_var_freq
    Float p_value
    File? roi_bed
  }

  call varscanSomatic {
    input:
    reference=reference,
    reference_fai=reference_fai,
    reference_dict=reference_dict,
    tumor_bam=tumor_bam,
    tumor_bam_bai=tumor_bam_bai,
    normal_bam=normal_bam,
    normal_bam_bai=normal_bam_bai,
    strand_filter=strand_filter,
    min_coverage=min_coverage,
    varscan_min_var_freq=varscan_min_var_freq,
    p_value=p_value,
    roi_bed=roi_bed
  }
}
