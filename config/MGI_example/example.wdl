version 1.0

workflow three_step {
  input {
    File procs
    String pattern
  }

  call ps
  call cgrep {
    input: in_file=ps.procs
  }
  call wc {
    input: in_file=ps.procs
  }
}

task ps {
  command {
    ps
  }
  
  runtime {
          docker_image: "ubuntu:xenial"
          cpu: "1"
          memory_gb: "4"
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
  }

  command {
    grep '${pattern}' ${in_file} | wc -l > output.result
  }
  
  output {
    File result = "output.result"
  }

  runtime {
          docker_image: "ubuntu:xenial"
          cpu: "1"
          memory_gb: "4"
          queue: "research-hpc"
          resource: "rusage[gtmp=10, mem=4000]"
          job_group: '/myjobgroup/'
  }
 
  
}
 
task wc {
  input {
    File in_file
  }

  command {
    cat ${in_file} | wc -l > output2.result
  }

  output {
    File result = "output2.result"
  }
  
  runtime {
          docker_image: "ubuntu:xenial"
          cpu: "1"
          memory_gb: "4"
          queue: "research-hpc"
          resource: "rusage[gtmp=10, mem=4000]"
          job_group: '/myjobgroup/'
  }
 
  
}
