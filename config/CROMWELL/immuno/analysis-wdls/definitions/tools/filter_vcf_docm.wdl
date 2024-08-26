version 1.0

task filterVcfDocm {
  input {
    File docm_raw_variants
    File normal_bam
    File tumor_bam
    Boolean filter_docm_variants
  }

  Int space_needed_gb = 10 + round(size(docm_raw_variants, "GB")*2 + size([normal_bam, tumor_bam], "GB"))
  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/samtools-cwl:1.16.1"
    memory: "4GB"
    disks: "local-disk ~{space_needed_gb} HDD"
  }

  String output_vcf_file = "docm_filtered_variants.vcf"
  String set_filter_flag = if filter_docm_variants then "1" else "0"
  command <<<
    /usr/bin/perl -e '
    use strict;
    use warnings;
    use feature qw(say);

    my $samtools = q(/opt/samtools/bin/samtools);
    my $normal_header_str = `$samtools view -H "~{normal_bam}" | grep "^\@RG" | head -n 1`;
    my $tumor_header_str  = `$samtools view -H "~{tumor_bam}" | grep "^\@RG" | head -n 1`;

    my ($normal_name) = $normal_header_str =~ /\@RG.+\tSM:([ -~]+)/;
    my ($tumor_name)  = $tumor_header_str =~ /\@RG.+\tSM:([ -~]+)/;

    unless ($normal_name and $tumor_name) {
      die "Failed to get normal_name: $normal_name from ~{normal_bam} AND tumor_name: $tumor_name from ~{tumor_bam}";
    }

    my $docm_vcf_fh;
    if("~{docm_raw_variants}" =~ /.gz$/) {
      open($docm_vcf_fh, "gunzip -c ~{docm_raw_variants} |") or die("could not open ~{docm_raw_variants} to read");
    } else {
      open($docm_vcf_fh, "~{docm_raw_variants}") or die("could not open ~{docm_raw_variants} to read");
    }
    open(my $docm_out_fh, ">", "~{output_vcf_file}") or die("could not open ~{output_vcf_file} for write");

    my ($normal_index, $tumor_index);

    while (<$docm_vcf_fh>) {
      chomp;
      if (/^##/) {
        say $docm_out_fh $_;
      }
      elsif (/^#CHROM/) {
        if (~{set_filter_flag}) {
          say $docm_out_fh q(##FILTER=<ID=DOCM_ONLY,Description="ignore Docm variants">);
        }
        my @columns = split /\t/, $_;
        my %index = (
        $columns[9]  => 9,
        $columns[10] => 10,
        );
        ($normal_index, $tumor_index) = map{$index{$_}}($normal_name, $tumor_name);
        unless ($normal_index and $tumor_index) {
          die "Failed to get normal_index: $normal_index for $normal_name AND tumor_index: $tumor_index for $tumor_name";
        }
        $columns[9]  = $normal_name;
        $columns[10] = $tumor_name;
        my $header = join "\t", @columns;
        say $docm_out_fh $header;
      }
      else {
        my @columns = split /\t/, $_;
        my @tumor_info = split /:/, $columns[$tumor_index];
        my ($AD, $DP) = ($tumor_info[1], $tumor_info[2]);
        next unless $AD;
        my @AD = split /,/, $AD;
        shift @AD; #the first one is ref count

        for my $ad (@AD) {
          if ($ad > 5 and $ad/$DP > 0.01) {
            my ($normal_col, $tumor_col) = map{$columns[$_]}($normal_index, $tumor_index);
            $columns[9]  = $normal_col;
            $columns[10] = $tumor_col;
            if (~{set_filter_flag}) {
              $columns[6] = q(DOCM_ONLY);
            }
            else {
              $columns[6] = q(.);
            }
            my $new_line = join "\t", @columns;
            say $docm_out_fh $new_line;
            last;
          }
        }
      }
    }

    close($docm_vcf_fh);
    close($docm_out_fh);
    '
    >>>

    output {
      File docm_filtered_variants = output_vcf_file
    }
}

workflow wf {
  input {
    File docm_raw_variants
    File normal_bam
    File tumor_bam
    Boolean filter_docm_variants
  }

  call filterVcfDocm {
    input:
    docm_raw_variants=docm_raw_variants,
    normal_bam=normal_bam,
    tumor_bam=tumor_bam,
    filter_docm_variants=filter_docm_variants
  }
}
