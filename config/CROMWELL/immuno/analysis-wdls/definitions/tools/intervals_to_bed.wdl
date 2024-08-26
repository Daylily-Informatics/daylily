version 1.0

task intervalsToBed {
  input {
    File interval_list
  }

  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "ubuntu:bionic"
    memory: "4GB"
  }

  command <<<
    /usr/bin/perl -e '
    use feature qw(say);

    for my $line (<>) {
      chomp $line;

      next if substr($line,0,1) eq q(@); #skip header lines

      my ($chrom, $start, $stop) = split(/\t/, $line);
      say(join("\t", $chrom, $start-1, $stop));
    }
    ' ~{interval_list} > interval_list.bed
  >>>

  output {
    File interval_bed = "interval_list.bed"
  }
}

workflow wf {
  input { File interval_list }
  call intervalsToBed {
    input: interval_list=interval_list
  }
}
