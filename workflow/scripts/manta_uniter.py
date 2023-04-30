"""
A script that will combine all of the candidateSV.vcf.gz with diploidSV.vcf.gz.
Prefrentially takes the diploidSV representation of alleles for QUAL/GT information.
"""

import sys
import gzip


def read(fn):
    """
    Returns {varkey: vcfline}
    """
    ret = {}
    sys.stderr.write("reading %s\n" % fn)
    cnt = 0
    with gzip.open(fn) as fh:
        for line in fh:
            line = line.decode()
            if line.startswith("#"):
                continue
            cnt += 1
            data = line.strip().split("\t")
            key = "%s:%s.%s.%s.%s" % (data[0], data[1], data[2], data[3], data[4])
            ret[key] = line
    sys.stderr.write("%d entries\n" % cnt)
    return ret


if __name__ == "__main__":
    # arg1 = diploidSV
    dip = read(sys.argv[1])
    # arg1 = candidateSV
    can = read(sys.argv[2])
    sys.stderr.write("header\n")
    # apparently need some of the other tags, too
    hlines = []
    with gzip.open(sys.argv[2]) as fh:
        for line in fh:
            line = line.decode()
            if line.startswith("##"):
                for i in [
                    "UPSTREAM_PAIR_COUNT",
                    "DOWNSTREAM_PAIR_COUNT",
                    "PAIR_COUNT",
                    "BND_PAIR_COUNT",
                ]:
                    if i in line:
                        hlines.append(line)
            else:
                break
    with gzip.open(sys.argv[1]) as fh:
        # get header
        for line in fh:
            line = line.decode()
            if line.startswith("##"):
                sys.stdout.write(line)
            else:
                sys.stdout.write("".join(hlines))
                sys.stdout.write(line)
                break
    sys.stderr.write("uniting\n")
    for key in can:
        if key in dip:
            sys.stdout.write(dip[key])
        else:
            data = can[key].strip().split("\t")
            data[5] = "0"
            data[6] = "PASS"
            data.append("GT\t./.\n")
            sys.stdout.write("\t".join(data))
    sys.stderr.write("finished\n")
