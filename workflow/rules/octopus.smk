import sys
import os

##### OCTOPUS- OUR snv CALLER
# ---------------------------
#
# One of the most important parts of the pipeline
# also one of the most complicated.
# Great docs can be found here: https://luntergroup.github.io/octopus/
# paper: https://www.nature.com/articles/s41587-021-00861-3

config["snv_pos_samps"] = {}

THREAD_PARTITION_MAP= {
'8' : 'i8,i64,i96',
'16'  : 'i8,i64,i96',
'32'  : 'i64,i96',
'64'  : 'i64,i96',
'96' : 'i96',
}


def get_oct_err_model(wildcards):

    # [library preparation]<.sequencer>. library preparation is selected from: PCR, PCR-FREE, or 10X. sequencer is selected from: HISEQ-2000, HISEQ-2500, HISEQ-4000, X10, NOVASEQ, BGISEQ-5000. 
    em=config["octopus"]["err_model"]  # ignoring this for now
    inst="NOVASEQ"  # See oct docs, but HISEQ, MISEQ
    lib_prep="PCR-FREE" # PCR or PCR-FREE
    if wildcards.sample in config['sample_info']:
        if 'instrument' in config['sample_info'][wildcards.sample]:
            inst=config['sample_info'][wildcards.sample]['instrument']
        if 'lib_prep' in config['sample_info'][wildcards.sample]:
            lib_prep = config['sample_info'][wildcards.sample]['lib_prep']      
    return f"{lib_prep}.{inst}"


def get_ploidy(wildcards):

    ploidy_str = " "
    try:
        if config["sample_info"][wildcards.sample]["biological_sex"].lower() in [
            "male",
            "m",
        ]:
            ploidy_str = ploidy_str + " --contig-ploidies X=1 Y=1 "
        elif config["sample_info"][wildcards.sample]["biological_sex"].lower() in [
            "female",
            "f",
        ]:
            ploidy_str = ploidy_str + " --contig-ploidies X=2 Y=0 "
        else:
            ploidy_str = ploidy_str + " --contig-ploidies X=2 Y=1 "
    except Exception as e:
        print(e, file=sys.stderr)

    return " "  # ploidy_str


# swap with conda directive to run
#    container:
#        "docker://dancooke/octopus"
#    container:
#        "docker://daylilyinformatics/octopus-skylake:0.7.4"


def get_ochrm_mod(wildcards):
    pchr=""
    if config['genome_build'] not in ['b37']:
        pchr="chr"
    ret_str = ""
    sl = wildcards.ochrm.split("-")
    sl2 = wildcards.ochrm.split("~")


    #from IPython import embed
    #embed()
    #raise

    if len(sl2) == 2:
        ret_str = pchr + wildcards.ochrm
    elif len(sl) == 1:
        ret_str = pchr + sl[0]
    elif len(sl) == 2:
        start = int(sl[0])
        end = int(sl[1])
        while start <= end:
            ret_str = str(ret_str) + " " + pchr + str(start)
            start = start + 1
    else:
        raise Exception(
            "oct chunks can only be one contiguous range per chunk : ie: 1-4 with the non numerical chrms assigned 23=X, 24=Y, 25=MT"
        )
    ret_str = ret_mod_chrm(ret_str)

    return ret_str


chrm_slots = {
    "1": "96",
    "2": "96",
    "3": "96",
    "4": "64",
    "5": "64",
    "6": "64",
    "7": "64",
    "8": "32",
    "9": "32",
    "10": "96",
    "11": "32",
    "12": "32",
    "13": "32",
    "14": "32",
    "15": "64",
    "16": "32",
    "17": "32",
    "18": "32",
    "19": "32",
    "20": "32",
    "21": "32",
    "22": "32",
    "23": "64",
    "24": "32",
    "25": "32",
}


def calc_oct_threads(wildcards, input):
    act_ochrm = str(wildcards.ochrm).split("-")[0].split("~")[0].split(":")[0]
    slot_n = 96
    if act_ochrm not in chrm_slots.keys():
        print(
            f" WEIRD CHROM NOT MATCHING OCT DICT:  {wildcards.ochrm}", file=sys.stderr
        )
    else:
        slot_n = chrm_slots[act_ochrm]

    return int(slot_n)


def calc_oct_partition(wildcards, input):
    act_ochrm = str(wildcards.ochrm).split("-")[0].split("~")[0].split(":")[0]
    slot_n = 96
    if act_ochrm not in chrm_slots.keys():
        print(
            f" WEIRD CHROM NOT MATCHING OCT DICT:  {wildcards.ochrm}", file=sys.stderr
        )
    else:
        slot_n = chrm_slots[act_ochrm]

    return str(THREAD_PARTITION_MAP[slot_n])


def calc_oct_x(wildcards):
    act_ochrm = str(wildcards.ochrm).split("-")[0].split("~")[0].split(":")[0]
    ret_x = " -X 4G  "
    if act_ochrm not in chrm_slots.keys():
        pass

    elif int(chrm_slots[act_ochrm]) == 96:
        ret_x = " -X 9G "
    elif int(chrm_slots[act_ochrm]) == 64:
        ret_x = " -X 8G "
    elif int(chrm_slots[act_ochrm]) == 32:
        ret_x = " -X  6G "
    elif int(chrm_slots[act_ochrm]) == 16:
        ret_x = " -X  5G "
    return f"{ret_x}"


def calc_oct_b(wildcards):
    act_ochrm = str(wildcards.ochrm).split("-")[0].split("~")[0].split(":")[0]
    ret_b = " -B 60M  "
    if act_ochrm not in chrm_slots.keys():
        pass
    elif int(chrm_slots[act_ochrm]) == 96:
        ret_x = " -B 180M "
    elif int(chrm_slots[act_ochrm]) > 64:
        ret_b = " -B 140M "
    elif int(chrm_slots[act_ochrm]) == 32:
        ret_b = " -B 100M "
    elif  int(chrm_slots[act_ochrm]) == 16:
        ret_b = " -B 75M "
    return f"{ret_b}"


def get_max_open_rds(wildcards):
    act_ochrm = str(wildcards.ochrm).split("-")[0].split("~")[0].split(":")[0]
    mor = "2500"
    if act_ochrm not in chrm_slots.keys():
        pass
    elif int(chrm_slots[act_ochrm]) == 96:
        mor = "60000"
    elif int(chrm_slots[act_ochrm]) == 64:
        mor = "38000"
    elif int(chrm_slots[act_ochrm]) == 32:
        mor="11000"
    elif int(chrm_slots[act_ochrm]) == 16:
        mor="5000"
    return f"{mor}"


def ret_oct_clust_samp(wildcards):
    return wildcards.sample + "_" + wildcards.ochrm


rule octopus:
    """https://github.com/luntergroup/octopus"""
    input:
        b=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
        bai=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam.bai",
        d=MDIR + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/{sample}.ready",
    output:
        #oct_tmpd=directory(MDIR + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/oct_tmp"),
        vcf=MDIR
        + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/{sample}.{alnr}.oct.{ochrm}.snv.vcf",
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/oct/log/vcfs/{sample}.{alnr}.oct.{ochrm}.snv.log",
    threads: calc_oct_threads  # config['octopus']['threads']
    container:
        "docker://daylilyinformatics/octopus-skylake:0.7.4" #conda: config['octopus']['env_yaml']
    priority: 45
    benchmark:
        repeat(
            MDIR + "{sample}/benchmarks/{sample}.{alnr}.oct.{ochrm}.bench.tsv",
            0
            if "bench_repeat" not in config["octopus"]
            else config["octopus"]["bench_repeat"],
        )
    resources:
        vcpu=calc_oct_threads,
        attempt_n=lambda wildcards, attempt:  (attempt + 0),
        partition="i192",
        threads=calc_oct_threads
    params:
        ts=calc_oct_threads,
        mor=get_max_open_rds,
        cluster_sample=ret_sample, #
        cluster_slots=calc_oct_threads,  #config['octopus']['threads'],  # calc_oct_threads,
        ochrm_mod=get_ochrm_mod,
        max_idel_err="  "
        if "max_idel_err" not in config["octopus"]
        else config["octopus"]["max_idel_err"],
        max_haplotypes=" "
        if "max_haplotypes" not in config["octopus"]
        else config["octopus"]["max_haplotypes"],
        max_open_read_files=config["octopus"]["max_open_read_files_modifier"],
        X=calc_oct_x,
        B=calc_oct_b,
        anno=config["octopus"]["anno"],
        target_working_memory=" ",  #config['octopus']['target_working_memory'],
        addl_options=" " if "addl_options" not in config["octopus"] else config["octopus"]["addl_options"],
        dup_read_detection=" ", # --duplicate-read-detection-policy AGGRESSIVE "        if "dup_policy" not in config["octopus"]        else config["octopus"]["dup_policy"]
        min_for_qual="  --min-forest-quality 0 " , #if "min_forest_quality" not in config["octopus"]        else config["octopus"]["min_forest_quality"],
        limit_calling_regions_to=config["octopus"]["regions_to_call"],
        tm=config["octopus"]["testing_model"],
        em=get_oct_err_model,
        fm=config['supporting_files']['files']['octopus']['forest_model_a']['name'],
        huref=config["supporting_files"]["files"]["octopus"]["huref"]["name"],
        skr=config['supporting_files']['files']['ucsc']['build_gaps']['name'],
        brt=" --bad-regon-tolerance NORMAL " if "brt" not in config["octopus"] else config["octopus"]["brt"],  
        ld_pre=config['octopus']['ld_pre'],
        mdir=MDIR,
        ploidy=" ", # get_ploidy,
        refcall=" ",
        tgt_working_mem=" --target-working-mem 6G "
        if "tgt_working_memory" not in config["octopus"]
        else config["octopus"]["tgt_working_memory"],
        octo_threads=calc_oct_threads,  #config['octopus']['threads'],
        model_posterior=" ",  #--model-posterior ALL ",
    shell:
        """
        export BRTOL="NORMAL";  

        export variable_args=" {params.addl_options} ";

        midel="{params.max_idel_err}";
        mhap="{params.max_haplotypes}";

        BRTL=" --bad-region-tolerance $BRTOL ";
        echo MOP {params.max_open_read_files};

        threads_flag=' --threads  ';

        timestamp=$(date +%Y%m%d%H%M%S);
        export TMPDIR=/fsx/scratch/octo_tmp_$timestamp;
        mkdir -p $TMPDIR;
        export APPTAINER_HOME=$TMPDIR;
        trap "rm -rf $TMPDIR" EXIT;
        
        export oochrm_mod=$(echo '{params.ochrm_mod}' | sed 's/~/\:/g' | perl -pe 's/(^23| 23)/ X/g;' | perl -pe 's/(^24| 24)/ Y/g;' | perl -pe 's/(^25| 25)/ MT/g;');

        ocmd="octopus -T $oochrm_mod $threads_flag     --reference {params.huref}  --reads {input.b}  $midel  $mhap  --annotations {params.anno}   --sequence-error-model {params.em}  {params.min_for_qual}   {params.tm} --skip-regions-file {params.skr}  $BRTL  {params.tgt_working_mem} $variable_args --max-open-read-files {params.mor} ";
        
        echo $ocmd 1&>2 > ocmd.log;

        $ocmd > {output.vcf};
        """


rule oct_sort_index_chunk_vcf:
    input:
        vcf=MDIR
        + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/{sample}.{alnr}.oct.{ochrm}.snv.vcf",
    priority: 46
    output:
        vcfsort=MDIR
        + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/{sample}.{alnr}.oct.{ochrm}.snv.sort.vcf",
        vcfgz=MDIR
        + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/{sample}.{alnr}.oct.{ochrm}.snv.sort.vcf.gz",
        vcftbi=MDIR
        + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/{sample}.{alnr}.oct.{ochrm}.snv.sort.vcf.gz.tbi",
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/oct/vcfs/{ochrm}/log/{sample}.{alnr}.oct.{ochrm}.snv.sort.vcf.gz.log",
    resources:
        vcpu=8,
        threads=8,
        partition="i192"
    params:
        cluster_sample=ret_sample,
    threads: 8 #config["config"]["sort_index_oct_chunk_vcf"]['threads']
    shell:
        """
        (rm {output} 1>  /dev/null  2> /dev/null )  || echo rmfailed > {log};
        (bedtools sort -header -i {input.vcf} > {output.vcfsort} 2>> {log}) || exit 1233;
        (
        bgzip {output.vcfsort};
        
        touch {output.vcfsort};
        tabix -f -p vcf {output.vcfgz};
        touch {output.vcftbi};
        {latency_wait};
        ls {output}; ) > {log} 2>&1 ;
        
        {latency_wait};
        """


localrules:
    oct_concat_fofn,


rule oct_concat_fofn:
    input:
        chunk_tbi=sorted(
            expand(
                MDIR
                + "{{sample}}/align/{{alnr}}/snv/oct/vcfs/{ochm}/{{sample}}.{{alnr}}.oct.{ochm}.snv.sort.vcf.gz.tbi",
                ochm=OCTO_CHRMS,            ),            key=lambda x: float(                str(x.replace("~", ".").replace(":", "."))               .split("vcfs/")[1]                .split("/")[0]                .split("-")[0]            ),        ),
    # This expand pattern is neat.  the escaped {} remain acting as a snakemake wildcard and expect to be derived from the dag, while th ochrm wildcard is effectively being constrained by the values in the OCTO_CHRMS array;  So you produce 1 input array of files for every sample+ochrm parir, with one list string/file name per array.  The rule will only begin when all array members are produced. It's then sorted by first octochrm so they can be concatenated w/out another soort as all the chunks had been sorted already.
    priority: 44
    output:
        fin_fofn=MDIR
        + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.concat.vcf.gz.fofn",
        tmp_fofn=MDIR        + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.concat.vcf.gz.fofn.tmp",
    threads: 2
    resources:
        threads=2
    params:
        fn_stub="{sample}.{alnr}.oct.",
        cluster_sample=ret_sample,
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.oct.concat.fofn.bench.tsv"
    conda:
        "../envs/vanilla_v0.1.yaml"
    log:
        MDIR + "{sample}/align/{alnr}/snv/oct/log/{sample}.{alnr}.oct.cocncat.fofn.log",
    shell:
        """
        (rm {output} 1> /dev/null  2> /dev/null ) || echo rmFailOK >> {log} && ls ./ >> {log};
        ### export LD_LIBRARY_PATH=$PWD/resources/libs/;
        for i in {input.chunk_tbi}; do
            ii=$(echo $i | perl -pe 's/\.tbi$//g'; );
            echo $ii >> {output.tmp_fofn};
        done;
        (workflow/scripts/sort_concat_chrm_list.py {output.tmp_fofn} {wildcards.sample}.{wildcards.alnr}.oct. {output.fin_fofn}) || echo "Python Script Error? CODE __ $? __" >> {log} && ls -lt {output} >> {log};

        """


rule oct_concat_index_chunks:
    input:
        fofn=MDIR
        + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.concat.vcf.gz.fofn",
    output:
        vcf=touch(
            MDIR + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.sort.vcf"
        ),
        vcfgz=touch(
            MDIR + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.sort.vcf.gz"
        ),
        vcfgztbi=touch(
            MDIR
            + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.sort.vcf.gz.tbi"
        ),
    threads: 8
    resources:
        vcpu=8,
        threads=8,
        partition="i192"
    priority: 47
    params:
        huref=config["supporting_files"]["files"]["huref"]["fasta"]["name"],
        cluster_sample=ret_sample,
    resources:
        attempt_n=lambda wildcards, attempt:  (attempt + 0)
    benchmark:
        MDIR + "{sample}/benchmarks/{sample}.{alnr}.oct.merge.bench.tsv"
    conda:
        config["octopus"]["oct_gather_env"]
    log:
        MDIR
        + "{sample}/align/{alnr}/snv/oct/log/{sample}.{alnr}.oct.snv.merge.sort.gatherered.log",
    shell:
        """
        (rm {output} 1> /dev/null  2> /dev/null ) || echo rmFAIL;
        mkdir -p $(dirname {log});
        ### export LD_LIBRARY_PATH=$PWD/resources/libs/;
        bcftools concat -a -d all --threads {threads} -f {input.fofn}  -O v -o {output.vcf};
        bcftools view -O z -o {output.vcfgz} {output.vcf};
        bcftools index -f -t --threads {threads} -o {output.vcfgztbi} {output.vcfgz};
        stats_f=$(echo "{output.vcfgz}.bcf.stats");
        bcftools stats -F {params.huref}  {output.vcfgz} > $stats_f;
        {latency_wait}; > {log} """


localrules:
    clear_combined_octovcf,


rule clear_combined_octovcf:  # TARGET:  clear combined octo vcf so the chunks can be re-evaluated if needed.
    input:
        expand(
            MDIR + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.sort.vcf.gz",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    priority: 42
    shell:
        "(rm {input}*  1> /dev/null  2> /dev/null ) || echo 'file not found for deletion: {input}';"


localrules:
    produce_oct_vcf,


rule produce_oct_vcf:  # TARGET: just gen octo calls
    input:
        expand(
            MDIR
            + "{sample}/align/{alnr}/snv/oct/{sample}.{alnr}.oct.snv.sort.vcf.gz.tbi",
            sample=SSAMPS,
            alnr=ALIGNERS,
        ),
    output:
        "gatheredall.oct",
    priority: 48
    log:
        "gatheredall.oct.log",
    shell:
        """( touch {output} ;

        {latency_wait}; ls {output} ) >> {log} 2>&1;
        """


localrules:
    oct_prep_chunkdirs,


rule oct_prep_chunkdirs:
    input:
        b=MDIR + "{sample}/align/{alnr}/{sample}.{alnr}.mrkdup.sort.bam",
    output:
        expand(
            MDIR + "{{sample}}/align/{{alnr}}/snv/oct/vcfs/{ochrm}/{{sample}}.ready",
            ochrm=OCTO_CHRMS,
        ),
    log:
        MDIR + "{sample}/align/{alnr}/snv/oct/logs/{sample}.{alnr}.chunkdirs.log",
    shell:
        """
        ( echo {output}  ;
        mkdir -p $(dirname {output} );
        touch {output};
        ls {output}; ) > {log} 2>&1;
        """
 
