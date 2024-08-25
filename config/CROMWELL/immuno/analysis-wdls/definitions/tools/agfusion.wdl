version 1.0

task agfusion {
  input {
    File fusion_predictions
    File agfusion_database
    Boolean annotate_noncanonical = false
    String output_dir = "agfusion_results"
  }

  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/agfusion:1.3.11-ensembl-105"
    memory: "32GB"
    cpu: 4
  }

  command <<<
    /usr/local/bin/agfusion batch -a starfusion --middlestar \
    -f ~{fusion_predictions} -db ~{agfusion_database} \
    ~{true="--noncanonical" false="" annotate_noncanonical} \
    -o ~{output_dir}

    cd ~{output_dir} || exit 1
    zip -r "../agfusion_results.zip" . 1> /dev/null
    cd .. || exit 1
  >>>

  output {
    File annotated_fusion_predictions_zip = "agfusion_results.zip"
  }
}
