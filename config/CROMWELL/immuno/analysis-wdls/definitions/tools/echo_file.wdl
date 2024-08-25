version 1.0

task echoFile {
  input {}

  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "ubuntu:bionic"
  }

  command <<<
    echo "TEST" > outfile.txt
  >>>

  output {
    File out = "outfile.txt"
  }
}
