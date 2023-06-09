import vcfpy
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import class_weight

def build_classifier(dnascope_vcf_path, octopus_vcf_path, truth_vcf_path, output_vcf_path):
    # Load the DNAscope VCF file
    dnascope_vcf = vcfpy.Reader.from_path(dnascope_vcf_path)

    # Load the Octopus VCF file
    octopus_vcf = vcfpy.Reader.from_path(octopus_vcf_path)

    # Load the truth VCF file
    truth_vcf = vcfpy.Reader.from_path(truth_vcf_path)

    # Create a dictionary of true variants
    from IPython import embed

    true_variants = {}
    for record in truth_vcf:
 
        if len(record.ALT) == 1 and len(record.REF) == len(record.ALT[0].value):
            key = (record.CHROM, record.POS)
            true_variants[key] = (record.REF, record.ALT[0])

    # Collect the variant features and labels
    X = []
    y = []

    def get_variant_features(rec):
        return np.array([rec.QUAL, rec.INFO.get('DP'), rec.INFO.get('AD')])

    from IPython import embed
    for vcf in [dnascope_vcf, octopus_vcf]:
        for record in vcf:

            if len(record.ALT) == 1 and len(record.REF) == len(record.ALT[0].value):
                key = (record.CHROM, record.POS)
                if key in true_variants:
                    X.append(get_variant_features(record))
                    y.append(true_variants[key] == (record.REF, record.ALT[0].value))

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Calculate class weights to adjust for the imbalance in the dataset
    class_weights = class_weight.compute_class_weight('balanced', np.unique(y_train), y_train)

    # Train a random forest classifier with a weighted loss function
    clf = RandomForestClassifier(n_estimators=100, random_state=42, class_weight=class_weights)
    clf.fit(X_train, y_train)

    # Predict the labels for the test set
    y_pred = clf.predict(X_test)

    # Calculate the accuracy and F1-score
    accuracy = np.mean(y_pred == y_test)
    f1_score = 2 * (np.sum(y_pred & y_test)) / (np.sum(y_pred) + np.sum(y_test))

    # Filter the DNAscope VCF file
    dnascope_filtered = []
    for record in dnascope_vcf:
        if len(record.ALT) == 1 and len(record.REF) == len(record.ALT[0].value):
            key = (record.CHROM, record.POS)
            if key in true_variants and clf.predict(get_variant_features(record))[0]:
                dnascope_filtered.append(record)
    dnascope_writer = vcfpy.Writer.from_path(output_vcf_path + '_dnascope', dnascope_vcf.header)
    for record in dnascope_filtered:
        dnascope_writer.write_record(record)
    dnascope_writer.close()

import sys

dnas=sys.argv[1]
octo=sys.argv[2]
truth=sys.argv[3]
outv=sys.argv[4]
build_classifier(dnas, octo, truth, outv)
