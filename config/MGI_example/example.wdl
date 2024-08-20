version 1.0

workflow three_step {
  input {
    File procs
    String pattern
    String docker_image
    String partition = "i4-5"
    Int cpu
    Int memory
  }

  call ps {
    input: 
      docker_image=docker_image,
      partition=partition,
      cpu=cpu,
      memory=memory
  }
  call cgrep {
    input: 
      in_file=procs,
      pattern=pattern,
      docker_image=docker_image,
      partition=partition,
      cpu=cpu,
      memory=memory
  }
  call wc {
    input: 
      in_file=procs,
      docker_image=docker_image,
      partition=partition,
      cpu=cpu,
      memory=memory
  }
}

task ps {
  input {
    String partition
    Int cpu
    Int memory
    String docker_image
  }
  command {
    ps
  }
  
  runtime {
          docker: docker_image
          cpu: 1
          memory: 4
          project: 'RandD'
          all_partitions: 'i4-5'
  }
 
  output {
    File procs = stdout()
  }
}
 
task cgrep {
  input {
    String pattern
    File in_file
    String docker_image
    String partition
    Int cpu
    Int memory
  }

  command {
    sleep 60 && grep '${pattern}' ${in_file} | wc -l > output.result
  }
  
  output {
    File result = "output.result"
  }

  runtime {
          docker: docker_image
          cpu: 1
          memory: 4
          project: 'RandD'
          all_partitions: 'i4-5'
  }
 
  
}
 
task wc {
  input {
    File in_file
    String docker_image
    String partition
    Int cpu
    Int memory
  }

  command {
    sleep 60 && cat ${in_file} | wc -l > output2.result
  }

  output {
    File result = "output2.result"
  }
  
  runtime {
          docker:  docker_image
          cpu: 1
          memory: 4
          partition: 'i4-5'
          project: 'RandD'
          all_partitions: 'i4-5'
  }
 
  
}
