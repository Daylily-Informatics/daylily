version 1.0

task phlat {
  input {
    String phlat_name = "phlat" 
    File cram
    File cram_crai
    File reference
    File reference_fai
    Int nthreads = 8
    Int mem = 20
    String index_dir = "" # optional if indexes are within run.b38.sh default dir
   }

  Int space_needed_gb = 10 + round(5*size([cram, cram_crai, reference, reference_fai], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "~{mem}GB"
    cpu: nthreads
    docker: "mgibio/phlat:1.1_withindex"
    disks: "local-disk ~{space_needed_gb} HDD"
    bootDiskSizeGb: 3*space_needed_gb
  }

  command <<<
    /bin/bash /usr/bin/run.b38.sh --tag ~{phlat_name} --bam ~{cram} --ref-fasta ~{reference} --rs-dir . --index-dir ~{index_dir}
  >>>

  output {
    File phlat_summary = phlat_name + "_HLA.sum"
  }
}

workflow wf {
  input {
    String? phlat_name
    File cram
    File cram_crai
    File reference
    File reference_fai
    Int? nthreads
    Int? mem
    String? index_dir 
  }
  call phlat {
    input:
    phlat_name=phlat_name,
    cram=cram,
    cram_crai=cram_crai,
    reference=reference,
    reference_fai=reference_fai,
    nthreads=nthreads,
    mem=mem,
    index_dir=index_dir,
  }
} 
