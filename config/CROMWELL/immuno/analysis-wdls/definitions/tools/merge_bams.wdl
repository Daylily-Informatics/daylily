version 1.0

task mergeBams {
  input {
    Array[File] bams
    Boolean sorted = false
    String name = "merged"
  }

  Int cores = 4
  Int space_needed_gb = 10 + round(4*size(bams, "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/bam-merge:0.1"
    memory: "8GB"
    cpu: cores
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outname = name + ".bam"
  command <<<
    perl <<'CODE'
    #!/usr/bin/perl
    use strict;
    use warnings;

    use Getopt::Std;
    use File::Copy;

    my $nthreads = ~{cores};
    my $outfilename = "~{outname}";
    my $sorted = ~{true=1 false=0 sorted};

    my @bams = ("~{sep="\", \"" bams}");
    die 'missing input bams' unless scalar(@bams);

    #if there is only one bam, just copy it and index it
    if (scalar(@bams) == 1) {
        copy($bams[0], $outfilename) or die 'failed to copy file:' . $!;
    } else {
        if ($sorted) {
            my $rv = system((qw(/usr/bin/sambamba merge -t)), $nthreads, $outfilename, @bams);
            $rv == 0 or die 'failed to merge with sambamba';
        } else { #unsorted bams, use picard
            my @args = (
                'OUTPUT=' . $outfilename,
                'ASSUME_SORTED=true',
                'USE_THREADING=true',
                'SORT_ORDER=unsorted',
                'VALIDATION_STRINGENCY=LENIENT',
                map { 'INPUT=' . $_ } @bams
            );
            my $rv = system((qw(java -jar -Xmx6g /opt/picard/picard.jar MergeSamFiles)), @args);
            $rv == 0 or die 'failed to merge with picard';
        }
    }
    if ($sorted) {
       my $rv = system((qw(/usr/bin/sambamba index)), $outfilename);
       $rv == 0 or die 'failed to index';
    }
    CODE
  >>>

  output {
    File merged_bam = outname
  }
}

workflow wf {
  input { Array[File] bams }
  call mergeBams { input: bams=bams }
}
