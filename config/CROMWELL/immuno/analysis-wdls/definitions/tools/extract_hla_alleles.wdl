version 1.0

task extractHlaAlleles {
  input {
    File optitype_file
    File phlat_file
  }

  Int space_needed_gb = 10 + round(size(optitype_file, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    memory: "2GB"
    docker: "ubuntu:xenial"
    disks: "local-disk ~{space_needed_gb} HDD"
  }
  
  # first, extract HLA class I from the optitype file
  # second, extract HLA class II from the phlat file
  # third, ensure there are only 2 fields of accuracy for alleles
 
  String outname = "hla_calls_newline.txt"
  String temp = "temp.txt"
  command <<<
    /usr/bin/awk '{FS="\t";getline;for(n=2;n<=NF-2;n++){if($n==""){}else{printf "HLA-"$n"\n"}}}' ~{optitype_file} > ~{temp}
    grep "HLA_D" ~{phlat_file} | /usr/bin/awk '{FS="\t";if($2==""){}else{printf $2"\n"};if($3==""){}else{printf $3"\n"}}' >> ~{temp}
    /usr/bin/awk -F":" '{print $1 ":" $2}' ~{temp} > ~{outname}
  >>>

  output {
    Array[String] allele_string = read_lines(outname)
    File allele_file = outname
  }
}

workflow wf {
  input {
    File optitype_file 
    File phlat_file 
  }
  call extractHlaAlleles { 
    input: 
    optitype_file=optitype_file,
    phlat_file=phlat_file 
  }
}
