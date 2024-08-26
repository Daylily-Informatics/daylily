version 1.0

task bisulfiteVcf2bed {
  input {
    File vcf
    File reference
    File reference_fai
    Boolean assay_non_cpg_sites
  }

  Int space_needed_gb = 10 + round(size([vcf, reference], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/biscuit:0.3.8"
    memory: "16GB"
    cpu: 2
    disks: "local-disk ~{space_needed_gb} HDD"
    bootDiskSizeGb: space_needed_gb
  }

  command <<<
    set -eou pipefail
    if [[ "~{assay_non_cpg_sites}" == "false" ]]
    then
        #Creates a gzipped bed and a bedgraph that leaves out MT, random, GL contigs, etc
        /usr/bin/biscuit vcf2bed -t cg -k 1 -e ~{vcf} | /usr/bin/biscuit mergecg ~{reference} /dev/stdin -k 2 |  tee >(/bin/gzip >cpgs.bed.gz) | cut -f 1-4 | sort -k1,1 -k2,2n -S 12G | /usr/bin/perl -ne 'print $_ if $_ =~ /^(chr)?[1-9]?[0-9|X|Y]\s/' >cpgs.bedgraph
    else
        /usr/bin/biscuit vcf2bed -t cg -k 1 -e ~{vcf} | /usr/bin/biscuit mergecg ~{reference} /dev/stdin -k 2 |  tee >(/bin/gzip >cpgs.bed.gz) | cut -f 1-4 | sort -k1,1 -k2,2n -S 12G | /usr/bin/perl -ne 'print $_ if $_ =~ /^(chr)?[1-9]?[0-9|X|Y]\s/' >cpgs.bedgraph
        /usr/bin/biscuit vcf2bed -t ch -k 1 -e ~{vcf} | awk '$6 == "CA"' | tee >(/bin/gzip >cpas.bed.gz) | cut -f 1-3,8 | sort -k1,1 -k2,2n -S 12G | /usr/bin/perl -ne 'print $_ if $_ =~ /^(chr)?[1-9]?[0-9|X|Y]\s/' >cpas.bedgraph
        /usr/bin/biscuit vcf2bed -t ch -k 1 -e ~{vcf} | awk '$6 == "CT"' | tee >(/bin/gzip >cpts.bed.gz) | cut -f 1-3,8 | sort -k1,1 -k2,2n -S 12G | /usr/bin/perl -ne 'print $_ if $_ =~ /^(chr)?[1-9]?[0-9|X|Y]\s/' >cpts.bedgraph
        /usr/bin/biscuit vcf2bed -t ch -k 1 -e ~{vcf} | awk '$6 == "CC"' | tee >(/bin/gzip >cpcs.bed.gz) | cut -f 1-3,8 | sort -k1,1 -k2,2n -S 12G | /usr/bin/perl -ne 'print $_ if $_ =~ /^(chr)?[1-9]?[0-9|X|Y]\s/' >cpcs.bedgraph

    fi
  >>>

  output {
    Array[File] methylation_bed = glob("*.bed.gz")
    Array[File] methylation_bedgraph = glob("*.bedgraph")
  }
}
