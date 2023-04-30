import vcf
import pandas as pd
import numpy as np

# Load the VCF file into a pandas dataframe
vcf_file = "HG002.GRCh38.final.vcf.gz"
vcf_reader = vcf.Reader(filename=vcf_file, compressed=True)
vcf_df = pd.DataFrame.from_records((r.INFO, r.FORMAT, *r) for r in vcf_reader)
vcf_df.columns = list(vcf_reader.infos.keys()) + list(vcf_reader.formats.keys()) + vcf_reader._column_headers

# Filter out low-quality variants
vcf_df = vcf_df[(vcf_df["QUAL"] >= 20) & (vcf_df["FILTER"] == "PASS")]

# Calculate correlation coefficients for selected INFO and FORMAT fields
corr_info = vcf_df[["INFO_field_1", "INFO_field_2", "INFO_field_3"]].corr()
corr_format = vcf_df[["FORMAT_field_1", "FORMAT_field_2", "FORMAT_field_3"]].corr()

# Identify correlations between fields and the "Truth_Sensitivity" field
corr_truth_info = vcf_df[["INFO_field_1", "INFO_field_2", "INFO_field_3", "Truth_Sensitivity"]].corr()["Truth_Sensitivity"][:-1]
corr_truth_format = vcf_df[["FORMAT_field_1", "FORMAT_field_2", "FORMAT_field_3", "Truth_Sensitivity"]].corr()["Truth_Sensitivity"][:-1]

# Identify fields with high correlations to the "Truth_Sensitivity" field
high_corr_fields_info = list(corr_truth_info[abs(corr_truth_info) >= 0.5].index)
high_corr_fields_format = list(corr_truth_format[abs(corr_truth_format) >= 0.5].index)

# Print the fields with high correlations to the "Truth_Sensitivity" field
print("INFO fields with high correlation to Truth_Sensitivity: ", high_corr_fields_info)
print("FORMAT fields with high correlation to Truth_Sensitivity: ", high_corr_fields_format)

# Filter out variants with low correlations to the "Truth_Sensitivity" field
vcf_df = vcf_df[(vcf_df["Truth_Sensitivity"] >= 0.5) | (vcf_df["Truth_Sensitivity"].isnull())]

# Write the filtered VCF file to a new file
vcf_writer = vcf.Writer(open("HG002.filtered.vcf", "w"), vcf_reader)
for record in vcf_df.to_dict("records"):
    info = {k: v for k, v in record.items() if k in vcf_reader.infos}
    fmt = {k: v for k, v in record.items() if k in vcf_reader.formats}
    vcf_writer.write_record(vcf.model._Call(record["CHROM"], record["POS"], record["ID"], record["REF"], record["ALT"], None, info, fmt))
vcf_writer.close()
