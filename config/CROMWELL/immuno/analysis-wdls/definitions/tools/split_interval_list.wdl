version 1.0

task splitIntervalList {
  input {
    File interval_list
    Int scatter_count
  }

  Int space_needed_gb = 10 + round(size(interval_list, "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "broadinstitute/picard:2.24.2"
    memory: "6GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  command <<<
    /usr/bin/perl -e '
    use File::Copy;
    use warnings;
    use strict;
    my $output_dir = $ENV{q(PWD)};
    my $interval_list = q(~{interval_list});
    my $scatter_count = ~{scatter_count};
    my $i = 1;
    if ($scatter_count == 1) {
        File::Copy::copy($interval_list,qq{$i.interval_list});
    } else {
        my $retval = system(q{/usr/bin/java}, q{-jar}, q{/usr/picard/picard.jar}, q{IntervalListTools}, q{OUTPUT=}.$output_dir, q{INPUT=}.$interval_list, q{SCATTER_COUNT=}. $scatter_count);
        exit $retval if $retval != 0;
        for (glob(qq{*/scattered.interval_list})) {
            #create unique names and relocate all the scattered intervals to a single directory
            File::Copy::move($_, qq{$i.interval_list});
            $i++
        }
    }
    '
  >>>

  output {
    Array[File] split_interval_lists = glob("*.interval_list")
  }
}

workflow wf {
  input {
    File interval_list
    Int scatter_count
  }

  call splitIntervalList {
    input:
    interval_list=interval_list,
    scatter_count=scatter_count
  }
}
