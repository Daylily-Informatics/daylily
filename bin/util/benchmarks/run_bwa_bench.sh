

# activate the env you wish and then souce this to run
# nc = [ none, amd, navops, taskset  ]
# stmem = 1000M 500M..  #M
# proc_sub=
# ref= navops, clvr
# mode: bwaonly, serial, combined
# streamer: pigz, zcat, cat
# proc_sub : on off
# K: 100000-10000000000

outfile_stub=$1

hyperfine -S bash \
	  --export-json=$outfile_stub.json \
	  --export-csv=$outfile_stub.csv \
	  --export-markdown=$outfile_stub.md \
	  --show-output \
	  -r 3 \
	  -L nc none \
	  -L bwamem bwa-mem2 \
	  -L bwa_threads 95 \
	  -L sort_threads 50 \
	  -L view_threads 10 \
	  -L stmem 1500M \
	  -L fq1  '/jmajor/PHYTO__RESOURCES/source_analyses_data/MGtestdata/NA224694/R0_EX4743_SQ00_0.R1.fastq.gz' \
	  -L fq2 '/jmajor/PHYTO__RESOURCES/source_analyses_data/MGtestdata/NA224694/R0_EX4743_SQ00_0.R2.fastq.gz' \
	  -L proc_sub on,serial \
	  -L mode combine \
	  -L ref navops \
	  -L streamer unpigz \
          -L K 60000000 \
	  "./benchmark_bwamem2a.sh {nc} {bwamem} {bwa_threads} {sort_threads} {view_threads} {stmem} {proc_sub} {fq1} {fq2} {mode} {ref} {streamer} {K}" 2>&1

# '/fsx/jmajor/avx_tests/AVX512/512img/main/.tiny_data/data/R0_EX4743_SQ0_0.R1.fastq.gz' \
# '/fsx/jmajor/avx_tests/AVX512/512img/main/.tiny_data/data/R0_EX4743_SQ0_0.R2.fastq.gz' \
# /jmajor/PHYTO__RESOURCES/source_analyses_data/MGtestdata/NA224694/R0_EX4743_SQ00_0.R2.fastq.gz' \                                                                
# /jmajor/PHYTO__RESOURCES/source_analyses_data/MGtestdata/NA224694/R0_EX4743_SQ00_0.R1.fastq.gz' \                                                                    
	  
echo done

