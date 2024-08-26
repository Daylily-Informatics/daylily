version 1.0

task cnvnator {
  input {
    File bam
    Int bin_size = 100
    Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
    File reference
    File reference_fai
    String sample_name
  }

  runtime {
    preemptible: 1
    maxRetries: 2
    docker: "mgibio/cnvnator-cwl:0.4"
    memory: "20GB"
    bootDiskSizeGb: 10
  }

  command <<<
    #set up the environment
    source /opt/root/bin/thisroot.sh
    set -eou pipefail
    # set vars
    BIN_SIZE="~{bin_size}"
    CHROMOSOMES="~{sep=" " chromosomes}"
    REFERENCE="~{reference}"
    SAMPLE="~{sample_name}"

    # create directory to store fasta files. CNVnator wants chromosomes in individual files
    # while with later versions of CNVnator these steps will accept a fasta.gz file the cnvnator2VCF.pl will not
    mkdir FASTA_CHRS
    awk 'BEGIN { CHROM="" } { if ($1~"^>") CHROM=substr($1,2); print $0 > "FASTA_CHRS/"CHROM".fa" }' "$REFERENCE"

    # extract read mapping from input bam(single sample)
    cnvnator -root "$SAMPLE.root" -tree "~{bam}" -chrom $CHROMOSOMES
    # generate read depth histogram
    cnvnator -root "$SAMPLE.root" -his "$BIN_SIZE" -d FASTA_CHRS/ -chrom $CHROMOSOMES
    # calculate statistics
    cnvnator -root "$SAMPLE.root" -stat "$BIN_SIZE" -chrom $CHROMOSOMES
    # read depth signal partitioning
    cnvnator -root "$SAMPLE.root" -partition "$BIN_SIZE" -chrom $CHROMOSOMES
    # cnv calling
    cnvnator -root "$SAMPLE.root" -call "$BIN_SIZE" -chrom $CHROMOSOMES > "$SAMPLE.CNVnator.cn"

    # convert to vcf
    cnvnator2VCF.pl -reference "$REFERENCE" "$SAMPLE.CNVnator.cn" FASTA_CHRS/ >  "$SAMPLE.CNVnator.vcf"
    exit 0
  >>>

  output {
    File vcf = sample_name + ".CNVnator.vcf"
    File root_file = sample_name + ".root"
    File cn_file = sample_name + ".CNVnator.cn"
  }
}
