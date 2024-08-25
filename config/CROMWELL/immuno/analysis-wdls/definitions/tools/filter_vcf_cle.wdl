version 1.0

task filterVcfCle {
  input {
    File vcf
    Boolean filter
  }

  Int space_needed_gb = 10 + round(size(vcf, "GB")*2)
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/samtools-cwl:1.16.1"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String do_filtering = if filter then "filter" else "0"
  command <<<
    /usr/bin/perl - "~{vcf}" "$PWD" "~{do_filtering}" <<'EOF'
#! /usr/bin/perl

#Copyright (C) 2015 Feiyu Du <fdu@genome.wustl.edu>
#              and Washington University The Genome Institute

#This script is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY or the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.


use strict;
use warnings;

use feature qw(say);

die "Provide vep annotated vcf and output_dir" unless @ARGV >= 2;
my ($annotated_vcf, $outdir, $filter_flag) = @ARGV;

my @coding_categories = qw(
    splice_acceptor_variant
    splice_donor_variant
    splice_region_variant
    transcript_ablation
    transcript_amplification
    stop_retained_variant
    stop_gained
    stop_lost
    start_lost
    frameshift_variant
    inframe_insertion
    inframe_deletion
    missense_variant
    protein_altering_variant
    initiator_codon_variant
    incomplete_terminal_codon_variant
    synonymous_variant
    coding_sequence_variant
);

open (my $annotated_vcf_fh, '-|', '/bin/gunzip', '-c', $annotated_vcf)
    or die "couldn't open $annotated_vcf to read";
open (my $annotated_filtered_vcf_fh, ">", "$outdir/annotated_filtered.vcf")
    or die "couldn't open annotated_filtered.vcf to write";

while (<$annotated_vcf_fh>) {
    chomp;
    if (/^#/) {
        say $annotated_filtered_vcf_fh $_;
    }
    else {
        if ($filter_flag and $filter_flag eq 'filter') {
            my @columns = split /\t/, $_;
            my ($ref, $alt, $info) = map{$columns[$_]}(3, 4, 7);
            my @alts = split /,/, $alt;
            my ($caller) = $info =~ /;set=(\S+?);/;

            if (length($ref) == 1 and length($alts[0]) == 1) { #snvs
                say $annotated_filtered_vcf_fh $_;
            }
            else {
                if ($caller =~ /docm/) {
                    say $annotated_filtered_vcf_fh $_;
                }
                else {
                    my ($csq) = $info =~ /;CSQ=(\S+)/;
                    if (grep {$csq =~ /$_/}@coding_categories) {
                        say $annotated_filtered_vcf_fh $_;
                    }
                }
            }
        }
        else {
            say $annotated_filtered_vcf_fh $_;
        }
    }
}

close $annotated_vcf_fh;
close $annotated_filtered_vcf_fh;
EOF

  >>>

  output {
    File cle_filtered_vcf = "annotated_filtered.vcf"
  }
}
