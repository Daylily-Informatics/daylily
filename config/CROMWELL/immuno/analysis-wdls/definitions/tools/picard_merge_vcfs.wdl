version 1.0

task picardMergeVcfs {
  input {
    String merged_vcf_basename = "merged"
    File? sequence_dictionary
    Array[File] vcfs
  }

  Int space_needed_gb = 10 + round(size(sequence_dictionary, "GB") + size(vcfs, "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "40GB"
    docker: "broadinstitute/gatk:4.1.8.1"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /usr/bin/java -Xmx38g -jar /gatk/gatk.jar MergeVcfs -O ~{merged_vcf_basename}.vcf.gz \
    ~{if defined(sequence_dictionary) then "-D ~{sequence_dictionary}" else ""} \
    ~{sep=" " prefix("-I ", vcfs)}
  >>>

  output {
    File merged_vcf = merged_vcf_basename + ".vcf.gz"
    File merged_vcf_tbi = merged_vcf_basename + ".vcf.gz.tbi"
  }
}


workflow wf {
  input {
    String merged_vcf_basename = "merged"
    File? sequence_dictionary
    Array[File] vcfs
  }

  call picardMergeVcfs {
    input:
    merged_vcf_basename=merged_vcf_basename,
    sequence_dictionary=sequence_dictionary,
    vcfs=vcfs
  }
}
