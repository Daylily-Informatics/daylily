#!/usr/bin/python
import os

for f in open('./fq_sstats.sort'):
         sl=f.rstrip().split('\t')
         fn = sl[0]
         if fn in ['file']:
                  continue
         nrd = sl[3]
         expected_cov = float(nrd)*float(150)/350000000.0

         scale_factor = 35.0/expected_cov
         fn2=fn.replace('.R1',',R2')
         if scale_factor > 1.0:
                  os.system(f'echo "THIS SAMPLE is below 35x {expected_cov}   {fn2}\n" >> ds_files/TOO_LOW_COVERAGE_SAMPLES.log')
         else:
                  x = fn.split('/')[2]
                  of = f"./ds_files/{x}.downsamp"
                  print(f"seqkit sample  -s 11  -j 10 -p {scale_factor} -o {of} {fn} ")
