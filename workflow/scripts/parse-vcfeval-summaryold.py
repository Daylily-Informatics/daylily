
""" Take rtg vcfeval summary.txt file, pull out the 'None' filtered data and regurigtate it with some add'l values calc'd"""
import os
import sys
import pandas as pd
import re
import math

cpus_div4 = os.cpu_count() // 4  # Integer division

summary_fh = open(sys.argv[1], "r")
sample = sys.argv[2]  # Sample name
tgt_region_bed = sys.argv[3]  # bed  file of the regions used in vcf eval
cmp_footprint = sys.argv[4]  # Replaced treatment
subset = None
alt_id = sys.argv[5]
new_summary_out_fh = open(sys.argv[6], "w")
allvar_mean_dp=None
try:
    allvar_mean_dp=int(float(str(sys.argv[7])))
except Exception as e:
    print(e,file=sys.stderr)
    allvar_mean_dp = -1

print(f"ARGS {sys.argv}", file=sys.stderr)
alnr=sys.argv[8]
snv_caller=sys.argv[9]



cov_bin=None
if allvar_mean_dp >=0 and allvar_mean_dp < 1:
    cov_bin=0
elif allvar_mean_dp > 0 and allvar_mean_dp < 3:
    cov_bin=1
elif allvar_mean_dp >= 3 and allvar_mean_dp < 5:
    cov_bin=3
elif allvar_mean_dp >= 5 and allvar_mean_dp < 7:
    cov_bin=5
elif allvar_mean_dp >= 7 and allvar_mean_dp < 10:
    cov_bin=7
elif allvar_mean_dp >= 10 and allvar_mean_dp < 15:
    cov_bin=10
elif allvar_mean_dp >= 15 and allvar_mean_dp < 20:
    cov_bin=15
elif allvar_mean_dp >= 20 and allvar_mean_dp < 25:
    cov_bin=20
elif allvar_mean_dp >=25 and allvar_mean_dp < 30:
    cov_bin=25
elif allvar_mean_dp >=30 and allvar_mean_dp < 35:
    cov_bin=30
elif allvar_mean_dp >=35 and allvar_mean_dp < 40:
    cov_bin=35
elif allvar_mean_dp >=40 and allvar_mean_dp < 45:
    cov_bin=40
elif allvar_mean_dp >=45 and allvar_mean_dp < 50:
    cov_bin=45
elif allvar_mean_dp >= 50:
    cov_bin=50
else:
    cov_bin=-2

print(f"VARCOV: {allvar_mean_dp} -- COVBINIS: {cov_bin}", file=sys.stderr)

even_newer_summary = "{0}/snv_{2}_{1}_concordance.mqc.tsv".format(
    os.path.dirname(sys.argv[6]), cmp_footprint, sample
)

cmd = "cat {0} | awk -F'\t' 'BEGIN{{SUM=0}}{{ SUM+=$3-$2 }}END{{print SUM}}'".format(
    tgt_region_bed
)
tgt_region_size = float(os.popen(cmd).readline())


new_summary_out_fh.write(
    "Sample\tCmpFootprint\tTN\tTP\tFP\tFN\tPrecision\tSensitivity-Recall\tSpecificity\tFDR\tFscore\tTgtRegionSize\tSubset\tAltID\tAllVarMeanDP\n"
)
ctr = 0
for i in summary_fh:
    l = i.lstrip(" ").rstrip()
    if ctr == 3:
        sl = re.split(r" +", l)
        threshold = sl[0]
        tp_baseline = float(sl[1])
        tp = float(sl[2])
        fp = float(sl[3])
        fn = float(sl[4])
        f_measure = float(sl[7])
        tn = float(tgt_region_size - (tp + fn))

        new_summary_out_fh.write(
            "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}".format(
                sample,
                cmp_footprint,
                tn,
                tp,
                fp,
                fn,
                tp / (tp + fp),
                tp / (tp + fn),
                tn / (tn + fp),
                fp / (tp + fp),
                f_measure,
                tgt_region_size,
                subset,
                alt_id,
                allvar_mean_dp
            )
        )
    ctr = 1 + ctr

def _proc_vcf(vcf_n):
    new_vcf_n = f"{vcf_n.replace('.','_')}_stripped.vcf.gz"
    ccmd = f" bcftools  annotate  -x INFO/datasets,INFO/platforms,INFO/callsetnames,INFO/platformnames,INFO/callable,INFO/datasetnames  --threads {cpus_div4} -O z -o {new_vcf_n} {vcf_n}; tabix -f {new_vcf_n}; "
    print(f"{ccmd}", file=sys.stderr)
    os.system(ccmd)
    return(new_vcf_n)

new_summary_out_fh.close()

tp_vcf = sys.argv[1].replace("summary.txt", "tp.vcf.gz")
tp_count = sys.argv[1].replace("summary.txt", "tp.count")
tp_cmd = "env python workflow/scripts/classify_var_by_type_size.py {0} {1} {2} {3} {4} {5}".format(
    _proc_vcf(tp_vcf), "TP", tgt_region_size, tp_count, sample, alt_id
)
print(tp_cmd, file=sys.stderr)
os.system(tp_cmd)

fp_vcf = sys.argv[1].replace("summary.txt", "fp.vcf.gz")
fp_count = sys.argv[1].replace("summary.txt", "fp.count")
fp_cmd = "env python workflow/scripts/classify_var_by_type_size.py {0} {1} {2} {3} {4} {5}".format(
    _proc_vcf(fp_vcf), "FP", tgt_region_size, fp_count, sample, alt_id
)
print(fp_cmd, file=sys.stderr)
os.system(fp_cmd)


fn_vcf = sys.argv[1].replace("summary.txt", "fn.vcf.gz")
fn_count = sys.argv[1].replace("summary.txt", "fn.count")
fn_cmd = "env python workflow/scripts/classify_var_by_type_size.py {0} {1} {2} {3} {4} {5}".format(
    _proc_vcf(fn_vcf), "FN", tgt_region_size, fn_count, sample, alt_id
)
print(fn_cmd, file=sys.stderr)
os.system(fn_cmd)

### Gather per-var type metrics, generate a more detailed view of the new_summary.txt
ds_var = {"TP": {}, "FP": {}, "FN": {}}

for cnt in [fn_count, fp_count, tp_count]:
    if os.path.exists(cnt):
        pass
    else:
        continue
    cnf_fh = open(cnt, "r")
    ctr = 0
    for c in cnf_fh:
        print(c, file=sys.stderr)
        c_sl = c.rstrip().split("\t")

        if ctr > 0:
            print(
                f"BBBB:{c_sl[0]}",file=sys.stderr
            )
            call_class = c_sl[1]
            print(f"XXXXX {c_sl} {c} {cnt}", file=sys.stderr)
            ds_var[call_class]["SNPts"] = int(c_sl[2])
            ds_var[call_class]["SNPtv"] = int(c_sl[3])
            ds_var[call_class]["INS_50"] = int(c_sl[4])
            ds_var[call_class]["INS_gt50"] = int(c_sl[5])
            ds_var[call_class]["DEL_50"] = int(c_sl[6])
            ds_var[call_class]["DEL_gt50"] = int(c_sl[7])
            ds_var[call_class]["Indel_50"] = int(c_sl[8])
            ds_var[call_class]["Indel_gt50"] = int(c_sl[9])
            ds_var[call_class]["tgtRegionSize"] = float(c_sl[10])

        elif ctr == 0:
            print("AAA", file=sys.stderr)

        else:
            print(c, file=sys.stderr)
            raise

        ctr += 1


print(ds_var, file=sys.stderr)

df = pd.DataFrame(ds_var)
try:
    df["TgtRegionSize"] = df["TP"]["tgtRegionSize"]
except Exception as e:
    df["TgtRegionSize"] = None
try:
    df = df.drop("tgtRegionSize")
except Exception as e:
    pass

df["TN"] = None
for i in df.iterrows():
    if math.isnan(i[1]["TP"]):
        i[1]["TP"] = 0.0
        df["TP"][i[0]] = 0.0
    if math.isnan(i[1]["FN"]):
        i[1]["FN"] = 0.0
        df["FN"][i[0]] = 0.0
    if math.isnan(i[1]["FP"]):
        i[1]["FP"] = 0.0
        df["FP"][i[0]] = 0.0

    df["TN"][i[0]] = float(i[1]["TgtRegionSize"]) - (
        float(i[1]["TP"]) + float(i[1]["FN"])
    )
df.loc["All"] = df.sum()

df["TgtRegionSize"]["All"] = tgt_region_size
df["TN"]["All"] = tgt_region_size - (df["TP"]["All"] + df["FN"]["All"])


df["Specificity"] = None
df["Sensitivity-Recall"] = None
df["FDR"] = None
df["PPV"] = None
df["Precision"] = None
df["Fscore"] = None

for i in df.iterrows():
    try:
        df["Precision"][i[0]] = float(i[1]["TP"]) / (
            float(i[1]["TP"]) + float(i[1]["FP"])
        )
    except:
        df["Precision"][i[0]] = None
    try:
        df["Sensitivity-Recall"][i[0]] = float(i[1]["TP"]) / (
            float(i[1]["TP"]) + float(i[1]["FN"])
        )
    except:
        df["Sensitivity-Recall"][i[0]] = None
    try:
        df["Specificity"][i[0]] = float(i[1]["TN"]) / (
            float(i[1]["TN"]) + float(i[1]["FP"])
        )
    except:
        df["Specificity"][i[0]] = None
    try:
        df["FDR"][i[0]] = float(i[1]["FP"]) / (float(i[1]["TP"]) + float(i[1]["FP"]))
    except:
        df["FDR"][i[0]] = None
    try:
        df["PPV"][i[0]] = 1.0 - float(
            float(i[1]["FP"]) / (float(i[1]["TP"]) + float(i[1]["FP"]))
        )
    except:
        df["PPV"][i[0]] = None
    try:  # TPTPFP = Precision TPTPFN--Sensitivity
        df["Fscore"][i[0]] = 2.0 * (
            (
                (float(i[1]["TP"]) / (float(i[1]["TP"]) + float(i[1]["FP"])))
                * (float(i[1]["TP"]) / (float(i[1]["TP"]) + float(i[1]["FN"])))
            )
            / (
                (float(i[1]["TP"]) / (float(i[1]["TP"]) + float(i[1]["FP"])))
                + (float(i[1]["TP"]) / (float(i[1]["TP"]) + float(i[1]["FN"])))
            )
        )
    except:
        df["Fscore"][i[0]] = None

df["mqc_id"] = f"{sample}-{alnr}-{snv_caller}-{subset}"
df["Sample"] = sample
df["AltId"] = alt_id
df["CmpFootprint"] = cmp_footprint
df["Subset"] = subset
df["AllVarMeanDP"] = allvar_mean_dp
df['CovBin'] = cov_bin
df['Aligner'] = alnr
df['SNVCaller'] = snv_caller
#print_cols = ['Sample'] + list(set(list(df.columns)) - set(['Sample']))

print_cols = ['mqc_id','Sample','TgtRegionSize','TN','FN','TP','FP','Fscore','Sensitivity-Recall','Specificity', 'FDR', 'PPV', 'Precision','AltId', 'CmpFootprint', 'AllVarMeanDP', 'CovBin', 'Aligner','SNVCaller']
df.to_csv(even_newer_summary, sep="\t", columns=print_cols)

os.system(f"perl -pi -e 's/^\t/SNPClass\t/g;' {even_newer_summary}")
