version 1.0

task freemix {
  input {
    File verify_bam_id_metrics
  }

  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "python:3.10"
  }

  command <<<
    python <<CODE
    with open("~{verify_bam_id_metrics}", "r") as f:
        header = f.readline().split("\t")
        if len(header) >= 7 and header[6] == "FREEMIX":
            with open("contamination.txt", "w") as out:
                out.write(f.readline().split("\t")[6])
    CODE
  >>>

  output {
    File? out = "contamination.txt"
  }
}

workflow wf {
  input {
    File verify_bam_id_metrics
  }

  call freemix {
    input:
    verify_bam_id_metrics=verify_bam_id_metrics
  }
}
