import vcfpy
import argparse

def calculate_concordance(truth_vcf_path, sample_vcf_path):
    truth_vcf = vcfpy.Reader.from_path(truth_vcf_path)
    sample_vcf = vcfpy.Reader.from_path(sample_vcf_path)

    truth_variants = {}
    for record in truth_vcf:
        if len(record.ALT) != 1:
            continue
        if len(record.REF) != len(record.ALT[0]):
            continue
        truth_variants[(record.CHROM, record.POS)] = (record.REF, record.ALT[0])

    true_positives = 0
    false_positives = 0
    false_negatives = 0
    total_snv_count = 0

    for record in sample_vcf:
        if len(record.ALT) != 1:
            continue
        if len(record.REF) != len(record.ALT[0]):
            continue
        variant_key = (record.CHROM, record.POS)
        if variant_key in truth_variants:
            total_snv_count += 1
            if truth_variants[variant_key] == (record.REF, record.ALT[0]):
                true_positives += 1
            else:
                false_positives += 1
                false_negatives += 1
        else:
            false_positives += 1

    sensitivity = true_positives / (true_positives + false_negatives)
    specificity = true_negatives / (true_negatives + false_positives)
    precision = true_positives / (true_positives + false_positives)
    fscore = 2 * (precision * sensitivity) / (precision + sensitivity)

    print(f'SNV concordance statistics:')
    print(f'--------------------------')
    print(f'Total SNVs: {total_snv_count}')
    print(f'True positives: {true_positives}')
    print(f'False positives: {false_positives}')
    print(f'False negatives: {false_negatives}')
    print(f'Sensitivity: {sensitivity}')
    print(f'Specificity: {specificity}')
    print(f'Precision: {precision}')
    print(f'F-score: {fscore}')

    snv_sizes = range(1, 51)
    for snv_size in snv_sizes:
        true_positives = 0
        false_positives = 0
        false_negatives = 0
        total_snv_count = 0

        for record in sample_vcf:
            if len(record.ALT) != 1:
                continue
            if len(record.REF) != snv_size or len(record.ALT[0]) != snv_size:
                continue
            variant_key = (record.CHROM, record.POS)
            if variant_key in truth_variants:
                total_snv_count += 1
                if truth_variants[variant_key] == (record.REF, record.ALT[0]):
                    true_positives += 1
                else:
                    false_positives += 1
                    false_negatives += 1
            else:
                false_positives += 1

        sensitivity = true_positives / (true_positives + false_negatives)
        specificity = true_negatives / (true_negatives + false
