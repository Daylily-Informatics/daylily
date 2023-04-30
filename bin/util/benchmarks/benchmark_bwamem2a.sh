$!/usr/bin/env bash



nc=$1
bwamem=$2
bwa_threads=$3
sort_threads=$4
view_threads=$5
stmem=$6
proc_sub=$7
fq1=$8
fq2=$9
mode=${10}
ref=${11}  
streamer=${12}   # cat, zcat , unpigz 
K=${13}

if [[ "$streamer" == "unpigz" ]]; then
    streamer="unpigz -c -q -- "
fi

ref_seq=xxxx
if [[ "$ref" == "navops" ]]; then
    ref_seq="/jmajor/PHYTO__RESOURCES/genomic_data/organism_refrences/human/human_g1k_v37_modified.fasta/bwa-mem2vanilla/human_g1k_v37_modified.fasta"
else
    ref_seq="/jmajor/PHYTO__RESOURCES/genomic_data/organism_refrences/human/human_g1k_v37_modified.fasta/bwa-mem2vanilla/human_g1k_v37_modified.fasta"
fi

numa=""
if  [[ "$nc" == "navops" ]]; then
    numa="numactl -m 0 -C 0-71,24-94 -- "
elif [[ "$NC" == "all" ]]; then
    numa="numactl --interleave=all -- "
elif [[ "$nc" == "amd"  ]]; then
    numa="numactl --m 0 -C 0-71,24-94 -- "
elif [[ "$nc" == "taskset" ]]; then
    numa="taskset -a -c 0-94 "
elif [[ "$nc" == "none" ]]; then
    numa=""
fi

in_fq_str=" $fq1 $fq2 "

if [[ "$proc_sub" == "on" ]]; then
    in_fq_str="<($streamer   $fq1 )   <($streamer   $fq2 )   "
else
    in_fq_str=" $fq1 $fq2 "
fi 

tdir="./a"

if [[ "$mode" == "combine" ]]; then
    a="$numa $bwamem mem   -R '@RG\\tID:12345\\tSM:{params.rgsm}\\PL:illumina\\tPU:FAKEc'  -a  -Y    -K $K    -k 19  -t $bwa_threads   $ref_seq     $in_fq_str         | samtools sort -l 0  -m $stmem -@  $sort_threads -T $tdir -O SAM  -         | samtools view -b -1  -h -@ $view_threads -O BAM  -o out.sort.bam - ; samtools index out.sort.bam ;"
    echo "$a"

    eval $a


elif [[ "$mode" == "bwaonly" ]]; then

    a="$numa $bwamem mem   -R '@RG\\tID:12345\\tSM:{params.rgsm}\\PL:illumina\\tPU:FAKEc'  -a  -Y  -o out.sam  -K $K    -k 19  -t $bwa_threads   $ref_seq     $in_fq_str"
    echo "$a"
    eval $a
    
elif [[ "$mode" == "serial" ]]; then  
    
    a="$numa   $bwamem mem  -R '@RG\\tID:1234FAKEf\\tSM:FAKEe\\tPL:FAKEd\\tPU:FAKEc\\tCN:FAKEb\\tPG:FAKEa'       -o out.sam      -a  -Y    -K $K   -k 19  -t $bwa_threads         $ref_seq          $in_fq_str ;    samtools sort -l 0  -m $stmem -@  $write_threads -T $tdir -O SAM -o ./out.sort.sam  out.sam  ;   samtools view -b -1  -h -@ $bwa_threads -O BAM  -o out.sort.bam  out.sort.sam ; samtools index out.sort.bam " 
    echo "$a"

    eval $a
    
fi

