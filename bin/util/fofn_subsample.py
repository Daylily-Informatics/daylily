#!/usr/bin/python
import os
import sys

ROOT="results/mod/b37/"

#.R1.fastq.gz
## File input in format  The sample field, if not endin in gz, will then try to find a file in: #results/mod/b37/<sample>/<sample.replace(_0,_1)+".R1.fastq.gz"  #! Where the replace(_0) is done so it only trims from the end

#subsamp_pct   sampleid
#1.0 RU0_EX4970-coriel_SQ6205_0
#0.825369543052459 RU0_EX4970-coriel_SQ6206_0
#0.761319032327168 RU0_EX4970-coriel_SQ6207_0

#results/mod/b37/RU0_EX4970-coriel_SQ6210_0/RU0_EX4970-coriel_SQ6210_1.R1.fastq.gz
for f in open(sys.argv[1]):

         sl=f.rstrip().split(',')

         subsamp_pct = sl[0]
         samp=sl[1]
         if subsamp_pct in ['subsamp_pct']:
                  continue
         else:
                  subsamp_pct = float(subsamp_pct)
         if subsamp_pct > 1.0:
                  print(f'Sample :: {samp} has a coverage < the subsample target and can not meet the downsample target, its sampling ratio is {subsamp_pct}', file=sys.stderr)

         sstub = "_".join(samp.split('_')[:-1])
         fnR1 = f"{ROOT}/{samp}/{sstub}_1.R1.fastq.gz"
         ofR1 = f"DOWNSAMPLED_OUTPUT/output/{sstub}_1.R1.downsmapled32.5.fastq"
         print(f"seqkit sample  -s 11  -j 77 -p {subsamp_pct} -o {ofR1} {fnR1} ")
         fnR2 = f"{ROOT}/{samp}/{sstub}_1.R2.fastq.gz"
         ofR2 = f"DOWNSAMPLED_OUTPUT/output/{sstub}_1.R2.downsmapled32.5.fastq"
         print(f"seqkit sample  -s 11  -j 77 -p {subsamp_pct} -o {ofR2} {fnR2}")


         ofR2
