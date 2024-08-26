version 1.0

task addStrelkaGt {
  input {
    File vcf
  }

  Int space_needed_gb = 10 + round(size(vcf, "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "ubuntu:bionic"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String outfile = "add_gt.vcf"
  command <<<
    /usr/bin/perl -e '
    use strict;
    use warnings;

    use feature qw(say);

    open(my $strelka_vcf_fh, "-|", "/bin/gunzip", "-c", "~{vcf}") or die("could not open ~{vcf} to read");
    open(my $add_gt_fh, ">", "~{outfile}") or die("could not open ~{outfile} for write");

    while (<$strelka_vcf_fh>) {
        chomp;
        if (/^##/) {
            say $add_gt_fh $_;
        }
        elsif (/^#/) { #COLUMN HEADER
            say $add_gt_fh q(##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">);
            say $add_gt_fh $_;
        }
        else {
            my @columns = split /\t/, $_;
            my ($ref, $alt, $info) = map{$columns[$_]}(3, 4, 7);
            my @alts = split /,/, $alt;
            my ($n_gt, $t_gt);
            if (length($ref) == 1 and length($alts[0]) == 1) {
                my ($n_gt_info, $n_gt_str, $t_gt_str) = $info =~ /NT=(\S+?);QSS.*SGT=(\S+?)\->(\S+?);/;
                unshift @alts, $ref;
                my %ids;
                my $id = 0;
                for my $base (@alts) {
                    $ids{$base} = $id;
                    $id++;
                }

                my $parsed_n_gt = parse_gt($n_gt_str, \%ids);
                if ($n_gt_info eq "ref") {
                    $n_gt = "0/0";
                } elsif (defined($parsed_n_gt)) {
                    $n_gt = $parsed_n_gt;
                } else {
                    my $id_keys = join(" ", sort keys(%ids));
                    print STDERR "parse_gt for n_gt=\"$n_gt_str\" and ids=\"$id_keys\" resulted in undefined\n";
                    next;
                }

                my $parsed_t_gt = parse_gt($t_gt_str, \%ids);
                if (defined($parsed_t_gt)) {
                    $t_gt = $parsed_t_gt;
                } else {
                    my $id_keys = join(" ", sort keys(%ids));
                    print STDERR "parse_gt for t_gt=\"$t_gt_str\" and ids=\"$id_keys\" resulted in undefined\n";
                    next;
                }
            }
            else {#INDEL
                my ($n_gt_info, $t_gt_info) = $info =~ /;NT=(\S+?);.*SGT.*\->(\S+?);/;
                my %gt_info = (
                    ref => q(0/0),
                    het => q(0/1),
                    hom => q(1/1),
                    conflict => q(./.),
                );
                $n_gt = $gt_info{$n_gt_info};
                $t_gt = $gt_info{$t_gt_info};
            }
            $columns[8]  = q(GT:).$columns[8];
            $columns[9]  = $n_gt . q(:) . $columns[9];
            $columns[10] = $t_gt . q(:) . $columns[10];
            my $new_str = join "\t", @columns;
            say $add_gt_fh $new_str;
        }
    }
    close($strelka_vcf_fh);
    close($add_gt_fh);
    sub parse_gt {
        my ($gt_str, $ids) = @_;
        my @gt_ids = map{$ids->{$_}}(split //, $gt_str);
        # Cannot parse if either side is undefined
        my $i;
        foreach $i (@gt_ids) {
            if (!defined($i)) {
                return undef;
            }
        }

        return join q(/), sort @gt_ids;
    }
    #SNV example
    #1       10231   .       C       A       .       QSS_ref NT=ref;QSS=1;QSS_NT=1;SGT=AC->AC;SOMATIC;TQSS=2;TQSS_NT=2       DP:FDP:SDP:SUBDP:AU:CU:GU:TU    32:4:8:0:0,3:28,60:0,0:0,1      84:6:69:0:7,21:71,192:0,0:0,1
    #INDEL example
    ##1     965051  .       ATGTGTG A       .       QSI_ref IC=5;IHP=2;NT=ref;QSI=1;QSI_NT=1;RC=8;RU=TG;SGT=het->het;SOMATIC;TQSI=1;TQSI_NT=1       DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50   8:8:6,6:0,0:2,4:10.3:0.00:0.00  18:18:8,8:5,6:5,8:21:0.25:0.00
    '
  >>>

  output {
    File processed_vcf = outfile
  }
}

workflow wf {
  input {
    File vcf
  }

  call addStrelkaGt {
    input: vcf=vcf
  }
}
