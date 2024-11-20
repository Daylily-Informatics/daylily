import os

# #### fastqc
# -----------
# github: https://github.com/s-andrews/FastQC


rule fastqc_subsampled:
    input:
        fqr1s=get_raw_R1s,
        fqr2s=get_raw_R2s,
    output:
        f"{MDIR}" + "{sample}/seqqc/fastqc/{sample}.fastqc.done",
    benchmark:
        f"{MDIR}" + "{sample}/benchmarks/{sample}.fastqc.bench.tsv"
    threads: config["fastqc"]["threads"]
    resources:
        vcpu=config["fastqc"]["threads"]
    params:
        tmp=f"{MDIR}" + "{sample}/seqqc/fastqc/tmp",
        tool_dir=f"{MDIR}" + "{sample}/seqqc/fastqc",
        cluster_sample=ret_sample,
        fastqc_threads=config["fastqc"]["fastqc_threads"],
        subsample_pct="0.25" if 'subsample_pct' not in config['fastqc'] else config['fastqc']['subsample_pct']
    log:
        f"{MDIR}" + "{sample}/logs/fastqc/{sample}.fastqc.log",
    conda:
        config["fastqc"]["env_yaml"]
    shell:
        """
        rm -rf {params.tool_dir} ;
        mkdir -p {params.tool_dir} ;
        mkdir -p {params.tmp}  ;
        #fastqc -o {params.tool_dir} -t {threads} -d {params.tmp}  <(seqkit sample --proportion {params.subsample_pct} <(seqfu interleave -1 <(unpigz -c -q -- {input.fqr1s}) -2 <(unpigz -c -q -- {input.fqr2s}) ) )  ;
        fastqc -o {params.tool_dir} -t {params.fastqc_threads} -d {params.tmp}  {input.fqr1s}  {input.fqr2s}
        touch {output};
        touch {output}_subsampled_at_{params.subsample_pct};
        {latency_wait}; ls {output};
        """
#
#
#
# Despite it's ubiquity, I have not found fastqc very useful over the years.  But, if any QC is to happen, fastqc
# must be in the mix as someone will ask for it, and many just assume it's needed for a complete assessment of
# sequence quality and be alarmed not to see it present.
# It is also not concerned with R1, R2 or PE reads....
# I have no problem throwing in more QC tools, except when they become expensive. And fastqc is slow to run.  Couple
# this with the fact it rarely (rarely) offers an actionable insight that other, faster, more complete tools would
# provide as well.... then it starts to feel like an anchor.
#
# It was s/w built for another decade, and is no suprise it's not very relevant feeling now-a-days.  This said, age of code, or recent release dates are not what I am judging here (fqastqc last released in 2020)-  Bwa has not had an update since 2017(!) and Heng Li even advocates for building tools that do not incline themselves to constant change (  His comments: https://lh3.github.io/2019/03/11/on-maintaining-bioinformatics-software """This is how I maintain my personal projects: I try not to maintain them. To achieve the goal, 1) I strive to simplify the user interface to reduce users' questions. 2) I make my projects independent of each other and of external libraries to avoid changes in dependencies failing my own projects. 3) If I perceive significant changes to an existing code base, I more often duplicate the code and start fresh (e.g. fermi vs fermi2 vs fermi-lite, and minimap vs minimap2). This way I can forget about compatibility and freely cut technical debts without affecting the stability of the previous version.)"""  -- in anycase, fastqc is more importantly not delivering enough to justify it's long runtimes (and by running longer, more likely to crash due to  background computing instabilities(bsi), which are present in even the most stable environments, it's a matter of frequency.
#
# Long story short, I would like to run the tool, get a meaningful output, and not spend 4 hrs doing so.  Solution, run fastqc on a subsampled set of the input fastqs.   Which is complicated by the fact I may have 1-many input files in the form of paired R1&R2 fq.gz's.  So, I need a tidy way to take  'n' pairs of fastq files and have fastqc just process 'm' percent of the reads, so that it runs more quickly and yields a useful report.  Further, the sampling should be randomly from all pairs of files- maintaining the mate pairing between files when sampled.

# A brute force solution might be:  interleave each pair into 1 fastq, subsample from the interleaved file to a new file, and then combine all of the new files to present the new combined file to fastq.
#
# This is a good solution, I did not like the large number of file manipulations it could involve.  I came up with this solution that acheives the desired result.
#

# The input functions reurn a sorted list of R1 file names and a second list of matched sorted R2 file names.
# seqfu interleave takes a -1 and -2 file and will interleae the two.  It will not accept  multiple files per argument.   but using <(unpigz -c -q -- [sorted list of R1s]) and <(unpigz -c -q [sorted list of R2s],  it is effectively presented with one file of all the list members concatenated.  seqfu interleave spits out a stream of interleaved fastq to seqkit sample, which subsamples to a specified % (i'm using 0.25 but this can be configurable) which itself it wrapped in <() to present to fastqc as a file to read


localrules:
    just_fastqc,


rule just_fastqc:
    input:
        expand(MDIR + "{sample}/seqqc/fastqc/{sample}.fastqc.done", sample=SAMPS),
    output:
        "fqc.done",
    shell:
        "touch {output}"
