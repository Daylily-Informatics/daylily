version 1.0

task optitypeDna {
  input {
    String optitype_name = "optitype"
    File cram
    File cram_crai
    File reference
    File reference_fai
    Int threads = 8
    Int mem = 50
  }

  Int space_needed_gb = 10 + round(5*size([cram, cram_crai, reference, reference_fai], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "~{mem}GB"
    cpu: threads 
    docker: "mgibio/immuno_tools-cwl:1.0.2"
    disks: "local-disk ~{space_needed_gb} HDD"
    bootDiskSizeGb: 3*space_needed_gb
  }

  command <<<
    /bin/bash /usr/bin/optitype_script_wdl.sh /tmp . \
    ~{optitype_name} ~{cram} ~{reference} ~{threads} ~{mem}
  >>>

  output {
    File optitype_tsv = optitype_name + "_result.tsv"
    File optitype_plot = optitype_name + "_coverage_plot.pdf"
  }
}

workflow wf {
  input {
    String? optitype_name
    File cram
    File cram_crai
    File reference
    File reference_fai
    Int? threads
    Int? mem
  }
  call optitypeDna {
    input:
    optitype_name=optitype_name,
    cram=cram,
    cram_crai=cram_crai,
    reference=reference,
    reference_fai=reference_fai,
    threads=threads,
    mem=mem
  }
}
