version 1.0

task filterSvVcfReadSupport {
  input {
    Float abundance_percentage = 0.1
    File input_vcf
    String output_vcf_name = "filtered_sv.vcf"
    Int? paired_count
    Int? split_count
    String vcf_source  # enum ["manta" "smoove"]
  }

  Int space_needed_gb = 10
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    docker: "bcftools-cwl:1.12"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -eou pipefail

    abundance=~{abundance_percentage}
    paired_count=~{paired_count}
    split_count=~{split_count}
    vcf_source="~{vcf_source}"

    if [ "$vcf_source" == "smoove" ]; then
        echo "Running filter for smoove vcf"
        filter_expression="((AS >= $split_count) && (AP >= $paired_count) && ((AP+AS) / (AP+RP+AS+RS) >= $abundance))"
    elif [ "$vcf_source" ==  "manta" ]; then
        echo "Running filter for manta vcf"
        filter_expression="((SR[0:*]=\".\" || (SR[0:1] >= $split_count)) && (PR[0:1] >= $paired_count) && ((SR[0:*]=\".\" && (PR[0:1] / (PR[0:0]+PR[0:1]) >= $abundance))|| ((SR[0:1]+PR[0:1]) / (SR[0:0]+SR[0:1]+PR[0:1]+PR[0:0]) >= $abundance)))"
    else
        echo "vcf source: '$vcf_source' is not supported for SV filtering"
        exit 1
    fi

    /opt/bcftools/bin/bcftools filter "~{input_vcf}" -o "~{output_vcf_name}" -i "$filter_expression"
  >>>

  output {
    File filtered_sv_vcf = output_vcf_name
  }
}
