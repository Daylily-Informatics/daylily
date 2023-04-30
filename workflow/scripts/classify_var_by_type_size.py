
""" 3yr+? OLD OLD SUP CODE PORTED OVER """
import os
import sys
import vcf

in_vcf = sys.argv[1]
call_class = sys.argv[2]
tgt_region_size = sys.argv[3]
count_out = sys.argv[4]
sample_mg = sys.argv[5]
sample_coriel = sys.argv[6]

in_vcf_sl = in_vcf.split(".")  # must be ./VCF path, not full
in_vcf_fh = open(in_vcf, "rb")

out_count_fh = open(count_out, "w")

snp_tv_vcf = "{0}.X_tv_snp.vcf".format(in_vcf)
snp_tv_count = 0
snp_ts_vcf = "{0}.X_ts_snp.vcf".format(in_vcf)
snp_ts_count = 0
del_50_vcf = "{0}.X_del_50.vcf".format(in_vcf)
del_50_count = 0
del_gt50_vcf = "{0}.X_del_gt50.vcf".format(in_vcf)
del_gt50_count = 0
ins_50_vcf = "{0}.X_ins_50.vcf".format(in_vcf)
ins_50_count = 0
ins_gt50_vcf = "{0}.X_ins_gt50.vcf".format(in_vcf)
ins_gt50_count = 0
indel_50_vcf = "{0}.X_indel_50.vcf".format(in_vcf)
indel_count_50 = 0
indel_gt50_vcf = "{0}.X_indel_gt50.vcf".format(in_vcf)

indel_count_gt50 = 0


vcf_reader = vcf.Reader(in_vcf_fh)

vcf_out_ds = {
    "classes": {
        "snp_ts": vcf.Writer(open(snp_ts_vcf, "w"), vcf_reader),
        "snp_tv": vcf.Writer(open(snp_tv_vcf, "w"), vcf_reader),
        "del_50": vcf.Writer(open(del_50_vcf, "w"), vcf_reader),
        "del_gt50": vcf.Writer(open(del_gt50_vcf, "w"), vcf_reader),
        "ins_50": vcf.Writer(open(ins_50_vcf, "w"), vcf_reader),
        "ins_gt50": vcf.Writer(open(ins_gt50_vcf, "w"), vcf_reader),
        "indel_50": vcf.Writer(open(indel_50_vcf, "w"), vcf_reader),
        "indel_gt50": vcf.Writer(open(indel_gt50_vcf, "w"), vcf_reader),
    }
}


ctr = 0
ctr2 = 0
for record in vcf_reader:
    try:
        vt = record.var_type
        vst = record.var_subtype

        if vt == "snp":
            if vst == "ts":
                vcf_out_ds["classes"]["snp_ts"].write_record(record)
                snp_ts_count += 1
            elif vst == "tv":
                vcf_out_ds["classes"]["snp_tv"].write_record(record)
                snp_tv_count += 1
        elif vt == "indel":
            if vst == "ins":
                if len(record.ALT) == 1:
                    if len(record.ALT[0]) < 51:
                        vcf_out_ds["classes"]["ins_50"].write_record(record)
                        ins_50_count += 1
                    elif len(record.ALT[0]) > 50 and len(record.ALT[0]) < 51:
                        vcf_out_ds["classes"]["ins_gt50"].write_record(record)
                        ins_gt50_count += 1
                    else:
                        vcf_out_ds["classes"]["ins_gt50"].write_record(record)
                        ins_gt50_count += 1
                else:
                    l = 0
                    for il in record.ALT:
                        if len(il) > l:
                            l = il
                    if l < 51:
                        vcf_out_ds["classes"]["ins_50"].write_record(record)
                        ins_50_count += 1
                    elif l > 50:
                        vcf_out_ds["classes"]["ins_gt50"].write_record(record)
                        ins_gt50_count += 1
                    else:
                        vcf_out_ds["classes"]["ins_gt50"].write_record(record)
                        ins_gt50_count += 1
            elif vst == "del":
                dellen = record.end - record.start
                if dellen < 51:
                    vcf_out_ds["classes"]["del_50"].write_record(record)
                    del_50_count += 1
                elif dellen > 50:
                    vcf_out_ds["classes"]["del_gt50"].write_record(record)
                    del_gt50_count += 1
                else:
                    vcf_out_ds["classes"]["del_gt50"].write_record(record)
                    del_gt50_count += 1
            else:
                # weird indel
                ilen = 0
                for a in record.ALT:
                    ilen = ilen + len(a) + (record.end - record.start)
                if ilen < 51:
                    vcf_out_ds["classes"]["indel_50"].write_record(record)
                    indel_count_50 += 1
                else:
                    vcf_out_ds["classes"]["indel_gt50"].write_record(record)
                    indel_count_gt50 += 1
        else:
            ilen = 0
            for a in record.ALT:
                ilen = ilen + len(a) + (record.end - record.start)
            if ilen < 51:
                vcf_out_ds["classes"]["indel_50"].write_record(record)
                indel_count_50 += 1
            else:
                vcf_out_ds["classes"]["indel_gt50"].write_record(record)
                indel_count_gt50 += 1
        ctr = ctr + 1
        if ctr > 10000:
            print(".", ctr2)
            ctr = 0
            ctr2 = ctr2 + 1
    except Exception as e:
        print("WHAT HAPPENED?")

        ilen = 0
        for a in record.ALT:
            ilen = ilen + len(a) + (record.end - record.start)
        print("!!!!!!! Error in parsing VCF  ", record, "--", e, " LEN: ", ilen)

        if ilen < 51:
            vcf_out_ds["classes"]["indel_50"].write_record(record)
            indel_count_50 += 1
        else:
            vcf_out_ds["classes"]["indel_gt50"].write_record(record)
            indel_count_gt50 += 1
            vcf_out_ds['classes']['indel'].write_record(record)


for i in vcf_out_ds["classes"]:
    vcf_out_ds["classes"][i].close()

os.system("bgzip -f {0}".format(snp_tv_vcf))
os.system("bgzip -f {0}".format(snp_ts_vcf))
os.system("bgzip -f {0}".format(indel_50_vcf))
os.system("bgzip -f {0}".format(indel_gt50_vcf))
os.system("bgzip -f {0}".format(ins_50_vcf))
os.system("bgzip -f {0}".format(ins_gt50_vcf))
os.system("bgzip -f {0}".format(del_50_vcf))
os.system("bgzip -f {0}".format(del_gt50_vcf))


os.system("tabix -f {0}.gz".format(snp_tv_vcf))
os.system("tabix -f {0}.gz".format(snp_ts_vcf))
os.system("tabix -f {0}.gz".format(indel_50_vcf))
os.system("tabix -f {0}.gz".format(indel_gt50_vcf))
os.system("tabix -f {0}.gz".format(ins_50_vcf))
os.system("tabix -f {0}.gz".format(ins_gt50_vcf))
os.system("tabix -f {0}.gz".format(del_50_vcf))
os.system("tabix -f {0}.gz".format(del_gt50_vcf))

out_count_fh.write(
    "Sample\tCallClass\tSNPts\tSNPtv\tIns50\tIns_gt50\tDel50\tDel_gt50\tIndel_50\tIndel_gt50\tTgtRegionSize\tAltName\n"
)

out_count_fh.write(
    "{10}\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{11}\n".format(
        call_class,
        snp_ts_count,
        snp_tv_count,
        ins_50_count,
        ins_gt50_count,
        del_50_count,
        del_gt50_count,
        indel_count_50,
        indel_count_gt50,
        tgt_region_size,
        sample_mg,
        sample_coriel,
    )
)

out_count_fh.close()
