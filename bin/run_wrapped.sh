bwa-mem2 mem  -Y -K 150000000000 -k 19 -t 254 -p \
 /fsx/data/genomic_data/organism_references/H_sapiens/b37/human_1kg_v37/bwa-mem2vanilla/human_g1k_v37.fasta \
   <(rclone cat awsdaylily:.bam | samtools view --input-fmt=BAM --output-fmt=BAM -u -h -@ 5 | bamtofastq  S=s O=o 02=o2   ) \
  | bamsort level=1 SO=coordinate blockmb=60000 index=1 indexfilename=testO.bam.bai \
    inputformat=sam outputformat=bam inputthreads=5 outputthreads=5 O=testO.bam \
    fixmates=1 calmdn=1 calmdnmreference=/fsx/data/genomic_data/organism_references/H_sapiens/b37/human_1kg_v37/bwa-mem2vanilla/human_g1k_v37.fasta \
     calmdnmrecompindetonly=1 adddupmarksupport=1 markduplicates=1 sortthreads=25
