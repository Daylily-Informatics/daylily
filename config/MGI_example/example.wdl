ask ps {
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
  String pattern
  File in_file
  
  command {
    grep '${pattern}' ${in_file} | wc -l
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
    Int count = read_int(stdout())
  }
  
}
 
task wc {
  File in_file
  
  command {
    cat ${in_file} | wc -l
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
    Int count = read_int(stdout())
  }
  
}
 
workflow three_step {
  call ps
  call cgrep {
    input: in_file=ps.procs
  }
  call wc {
    input: in_file=ps.procs
  }
}