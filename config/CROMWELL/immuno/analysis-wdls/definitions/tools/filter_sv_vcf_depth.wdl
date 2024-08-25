version 1.0

task filterSvVcfDepth {
  input {
    File input_vcf
    String output_vcf_name = "filtered_cnv.vcf"
    Float deletion_depth = 0.75
    Float duplication_depth = 1.25
    String vcf_source  # enum ["cnvnator" "cnvkit" "duphold"]
  }

  Int space_needed_gb = 10 + round(2*size(input_vcf, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "4GB"
    docker: "mgibiobcftools-cwl:1.12"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    set -eou pipefail

    DEL_DEPTH="~{deletion_depth}"
    DUP_DEPTH="~{duplication_depth}"
    VCF_SOURCE="~{vcf_source}"

    if [ "$VCF_SOURCE" == "cnvnator" ]; then
        echo "Running filter for cnvnator vcf"
        filter_expression="(natorRD>$DUP_DEPTH | natorRD<$DEL_DEPTH)"
    elif [ "$VCF_SOURCE" ==  "cnvkit" ]; then
        echo "Running filter for cnvkit vcf"
        filter_expression="(FOLD_CHANGE>$DUP_DEPTH | FOLD_CHANGE<$DEL_DEPTH)"
    elif [ "$VCF_SOURCE" ==  "duphold" ]; then
        echo "Running filter for vcf annotated by duphold"
        filter_expression="((SVTYPE='DEL' & FMT/DHFFC[0] < $DEL_DEPTH) | (SVTYPE='DUP' & FMT/DHBFC[0] > $DUP_DEPTH) | SVTYPE='BND' | SVTYPE='INS' | SVTYPE='INV')"
    else
        echo "vcf source: '$VCF_SOURCE' is not supported for SV filtering"
        exit 1
    fi
    /opt/bcftools/bin/bcftools filter "~{input_vcf}" -o "~{output_vcf_name}" -i "$filter_expression"
  >>>

  output {
    File filtered_sv_vcf = output_vcf_name
  }
}
