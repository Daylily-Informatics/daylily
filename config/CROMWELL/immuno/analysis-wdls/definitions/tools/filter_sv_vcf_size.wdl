version 1.0

task filterSvVcfSize {
  input {
    File input_vcf
    String output_vcf_name = "filtered_sv_size.vcf"
    Int sv_size = 50000
    String size_method  # enum ["max_len" "min_len"]
  }

  Int space_needed_gb = 10
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    docker: "mgibio/bcftools-cwl:1.12"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -eou pipefail

    SV_SIZE=~{sv_size}
    FILTER_METHOD="~{size_method}"

    if [ "$FILTER_METHOD" == "max_len" ]; then
        echo "Running filter for max size svs"
        filter_expression="ABS(SVLEN) <= $SV_SIZE"
    elif [ "$FILTER_METHOD" ==  "min_len" ]; then
        echo "Running filter for min size svs"
        filter_expression="ABS(SVLEN) >= $SV_SIZE"
    else
        echo "Filter method: '$FILTER_METHOD' is not supported for size SV filtering"
        exit 1
    fi
    /opt/bcftools/bin/bcftools filter "~{input_vcf}" -o "~{output_vcf_name}" -i "$filter_expression"
  >>>

  output {
    File filtered_sv_vcf = output_vcf_name
  }
}
