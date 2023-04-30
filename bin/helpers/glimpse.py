import numpy as np
import pandas as pd
from scipy.stats import chi2

# Define the GLIMPSE imputation function
def glimpse_impute(genotypes, haplotypes, reference_alleles, alternate_alleles):
    # Convert genotypes to allele dosages
    dosage = genotypes / 2.0

    # Compute the reference and alternate allele frequencies
    p_ref = np.sum(haplotypes * reference_alleles, axis=1) / (2 * haplotypes.shape[1])
    p_alt = 1 - p_ref

    # Compute the expected genotype frequencies
    exp_genotypes = np.outer(p_ref, p_ref) * 4 + np.outer(p_ref, p_alt) * 2 + np.outer(p_alt, p_alt)

    # Compute the variance of the expected genotype frequencies
    var_genotypes = exp_genotypes * (1 - exp_genotypes) * (2 * haplotypes.shape[1] + 1) / (2 * haplotypes.shape[1])

    # Compute the likelihood of each genotype
    genotype_likelihoods = np.zeros((genotypes.shape[0], 3))
    for i in range(genotypes.shape[0]):
        for j in range(3):
            x = dosage[i]
            mu = p_ref[i] * 2 + p_alt[i] * j
            var = var_genotypes[i]
            if var == 0:
                if x == mu:
                    genotype_likelihoods[i, j] = 1
                else:
                    genotype_likelihoods[i, j] = 0
            else:
                z = (x - mu) / np.sqrt(var)
                genotype_likelihoods[i, j] = np.exp(-z ** 2 / 2) / np.sqrt(2 * np.pi * var)

    # Normalize the genotype likelihoods
    normalizer = np.sum(genotype_likelihoods, axis=1)
    for i in range(genotypes.shape[0]):
        genotype_likelihoods[i, :] /= normalizer[i]

    # Compute the imputed genotypes
    imputed_genotypes = np.zeros_like(genotypes)
    for i in range(genotypes.shape[0]):
        for j in range(3):
            imputed_genotypes[i, j] = np.sum(haplotypes[:, j] * genotype_likelihoods[i, :])

    # Choose the imputed genotype with the highest likelihood
    imputed_genotypes = np.argmax(imputed_genotypes, axis=1)

    # Convert imputed genotypes back to alleles
    imputed_alleles = np.where(imputed_genotypes == 0, reference_alleles, alternate_alleles)
    imputed_alleles = np.where(imputed_genotypes == 2, alternate_alleles, imputed_alleles)

    return imputed_alleles

# Load the input data
genotypes = pd.read_csv('genotypes.csv').values
haplotypes = pd.read_csv('haplotypes.csv').values
reference_alleles = pd.read_csv('reference_alleles.csv').values.flatten()
alternate_alleles = pd.read_csv('alternate_alleles.csv').values.flatten()

# Impute the missing genotypes using GLIMPSE
imputed_alleles = glimpse_impute(genotypes, haplotypes, reference_alle
