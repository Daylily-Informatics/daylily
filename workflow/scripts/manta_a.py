import os
import sys

log = sys.argv[1]
vcf = sys.argv[2]
sortf = sys.argv[3]
work_dir = sys.argv[4]
vcfgz = sys.argv[5]
try:
    os.system(
        "python workflow/scripts/manta_uniter.py {0}/results/variants/diploidSV.vcf.gz {0}/results/variants/candidateSV.vcf.gz > {1}".format(
            work_dir, vcf
        )
    )
    os.system("bedtools sort -header -i {0} > {1} ".format(vcf, sortf, log))
    os.system(" bgzip {0} ".format(sortf, log, work_dir))
    os.system("tabix -f {0} ".format(vcfgz, log))
except Exception as e:
    os.system("rm {0} {1} {2}".format(vcf, sortf, vcfgz))
    print "MANNTA FAILED - was there enough coverage?", e
    raise
