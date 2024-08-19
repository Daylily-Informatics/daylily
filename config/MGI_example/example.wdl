version 1.0

workflow three_step {
  input {
    File procs
    String pattern
    String docker_image
    String partition
    Int cpu
    Int memory
  }

  call ps
  call cgrep {
    input: 
      in_file=procs,
      pattern=pattern
  }
  call wc {
    input: in_file=procs
  }
}

task ps {
  command {
    ps
  }
  
  runtime {
          docker: docker_image
          cpu: "1"
          memory: 4
          queue: "research-hpc"
          resource: "rusage[gtmp=10, mem=4000]"
          job_group: '/myjobgroup/'
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
          cpu: "1"
          memory: 4
          queue: "research-hpc"
          resource: "rusage[gtmp=10, mem=4000]"
          job_group: '/myjobgroup/'
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
          cpu: "1"
          memory: 4
          queue: "research-hpc"
          resource: "rusage[gtmp=10, mem=4000]"
          job_group: '/myjobgroup/'
  }
 
  
}
