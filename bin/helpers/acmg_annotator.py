
# conda create -n ACMG -c bioconda pyvcf

import sys
import argparse
import vcf as pyvcf
import requests

# Define the ACMG guidelines
ACMG_GUIDELINES = {
    'PVS1': {'criteria': 'null_variant'},
    'PS1': {'criteria': 'same_aa_change'},
    'PS2': {'criteria': 'de_novo_confirmed'},
    'PS3': {'criteria': 'well-established_function'},
    'PS4': {'criteria': 'prevalence'},
    'PM1': {'criteria': 'mut_hotspot'},
    'PM2': {'criteria': 'absent_in_controls'},
    'PM3': {'criteria': 'for_recessive'},
    'PM4': {'criteria': 'protein_length_change'},
    'PM5': {'criteria': 'novel_missense'},
    'PP1': {'criteria': 'co-segregated'},
    'PP2': {'criteria': 'missense_in_motif'},
    'PP3': {'criteria': 'multiple_computational'},
    'PP4': {'criteria': 'patient_phenotype'},
    'PP5': {'criteria': 'reputable_source'},
    'BA1': {'criteria': 'stand-alone'},
    'BS1': {'criteria': 'allele_data'},
    'BS2': {'criteria': 'observed_expected'},
    'BS3': {'criteria': 'in_vitro'},
    'BS4': {'criteria': 'non-segregation'},
    'BP1': {'criteria': 'missense_not_motif'},
    'BP2': {'criteria': 'observed_ctrl'},
    'BP3': {'criteria': 'non-pathogenic'},
    'BP4': {'criteria': 'other_database'},
    'BP5': {'criteria': 'functional_data'},
    'BP6': {'criteria': 'rebuttal'},
    'BP7': {'criteria': 'synonymous'}
}

def annotate_vcf(input_vcf, output_vcf):
    vcf_reader = pyvcf.Reader(open(input_vcf, 'r'))

    print(input_vcf, output_vcf)
    vcf_writer = pyvcf.Writer(open(output_vcf, 'w'), vcf_reader)

    # Iterate through each variant in the input VCF file
    for record in vcf_reader:
        info = record.INFO
        from IPython import embed
        embed()
        variant = info['ANN'][0]
        
        # Annotate the variant according to the ACMG guidelines
        for key, value in ACMG_GUIDELINES.items():
            criteria = value['criteria']
            
# Check if the variant meets the ACMG criteria and add it to the INFO field
            if criteria in variant:
                if 'ACMG' in info:
                    info['ACMG'].append(key)
                else:
                    info['ACMG'] = [key]
            
        # Write the annotated record to the output VCF file
        vcf_writer.write_record(record)

    # Close the output VCF file
    vcf_writer.close()

def main():
    # Set up command-line arguments
    parser = argparse.ArgumentParser(description='Annotate VCF with ACMG guidelines.')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file.')
    parser.add_argument('-o', '--output', required=True, help='Output VCF file.')

    # Parse command-line arguments
    args = parser.parse_args()

    # Annotate the input VCF file according to the ACMG guidelines
    annotate_vcf(args.input, args.output)

if __name__ == '__main__':
    main()


# This script takes an input human genome VCF file and annotates the variants according to the ACMG guidelines. To run the script, use the following command:
# bash
# python acmg_annotator.py -i input.vcf -o output.vcf

# Replace `input.vcf` with the path to your input VCF file, and `output.vcf` with the desired path for the output VCF file containing the annotated variants.
