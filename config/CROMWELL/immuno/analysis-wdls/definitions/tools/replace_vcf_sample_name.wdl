version 1.0

task replaceVcfSampleName {
  input {
    File input_vcf
    String sample_to_replace
    String new_sample_name
  }

  Int space_needed_gb = 10 + round(size(input_vcf, "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "8GB"
    docker: "mgibio/bcftools-cwl:1.12"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String basen = "renamed." + basename(input_vcf)
  command <<< #!/bin/bash
    set -eou pipefail
    #escape spaces, otherwise bcftools will try to use them as a delimiter
    #triple backslash to escape within backticks and then again within sed
    old_name=$(echo "~{sample_to_replace}" | sed 's/ /\\\ /g')
    new_name=$(echo "~{new_sample_name}" | sed 's/ /\\\ /g')

    echo "$old_name $new_name" > sample_update.txt
    /opt/bcftools/bin/bcftools reheader -s sample_update.txt -o "~{basen}" "~{input_vcf}"
  >>>

  output {
    File renamed_vcf = "renamed." + basename(input_vcf)
  }
}
