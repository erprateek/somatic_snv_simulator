import numpy as np
import argparse
import psutil

def get_available_memory():
    mem = psutil.virtual_memory()
    return mem.available

def simulate_somatic_mutations(num_pileups, depth, error_rate, allele_fraction, alt_threshold, fv):
    # Draw total read depth from a Poisson distribution
    total_depth            = np.random.poisson(depth, size=num_pileups)

    # Draw sequencing noise from a Binomial distribution
    noisy_alt_reads            = np.random.binomial(total_depth, error_rate)

    # Generate reads having mutation - Using 1 for True and 0 for False
    has_mutation           = np.random.choice([True, False], size=num_pileups, p=[fv, 1-fv])

    somatic_mutation_reads = np.random.binomial(total_depth, allele_fraction)*has_mutation

    total_alt_reads = noisy_alt_reads + somatic_mutation_reads

    return (total_alt_reads, has_mutation)


def calculate_confusion_matrix(labels, predictions):
    true_positive  = np.sum(labels & predictions)
    false_positive = np.sum(~labels & predictions)
    true_negative  = np.sum(~labels & ~predictions)
    false_negative = np.sum(labels & ~predictions)

    confusion_matrix = np.array([[true_negative, false_positive], [false_negative, true_positive]])
    return confusion_matrix

def calculate_ppa(confusion_matrix):
    true_positive = confusion_matrix[1, 1]
    false_negative = confusion_matrix[1, 0]

    ppa = true_positive / (true_positive + false_negative)
    return round(ppa*100, 2)

def calculate_ppv(confusion_matrix):
    true_positive = confusion_matrix[1, 1]
    false_positive = confusion_matrix[0, 1]

    ppv = true_positive / (true_positive + false_positive)
    return round(ppv*100,2)

def calculate_specificity(confusion_matrix):
    true_negative = confusion_matrix[0, 0]
    false_positive = confusion_matrix[0, 1]

    specificity = true_negative / (true_negative + false_positive)
    return round(specificity*100, 2)

def run_simulation_and_gather_metrics(num_pileups, depth, error_rate, allele_fraction, alt_threshold, fv):
    # Simulate somatic mutations
    pileups, labels = simulate_somatic_mutations(num_pileups,
                                                 depth,
                                                 error_rate,
                                                 allele_fraction,
                                                 alt_threshold,
                                                 fv)
    predictions = pileups >= alt_threshold
    confusion_matrix = calculate_confusion_matrix(labels, predictions)

    # Calculate PPA, PPV, and Specificity
    ppa = calculate_ppa(confusion_matrix)
    ppv = calculate_ppv(confusion_matrix)
    specificity = calculate_specificity(confusion_matrix)

    # Print the results
    return (confusion_matrix, ppa, ppv, specificity)

def main():
    parser = argparse.ArgumentParser("Monte-Carlo simulation of somatic mutations in NGS data.")
    parser.add_argument('-k','--num_pileups',          type=int,   help='Number of pileups to simulate')
    parser.add_argument('-D','--mean_total_depth',     type=float, help='Mean of the Poisson distribution, from which depth values are drawn')
    parser.add_argument('-e','--mean_seq_error_rate',  type=float, help='Average sequencing error rate. The number of ALT reads due to noise is drawn from BD(D,e)')
    parser.add_argument('-f','--frac_pileups_with_var',type=float, help='The fraction of k pileups with a somatic variant')
    parser.add_argument('-a','--af',                   type=float, help='The intrinsic Allele Fraction of somatic variants. The number of ALT reads due to a somatic variant is drawn from BD(D,AF)')
    parser.add_argument('-t','--alt_threshold',        type=int,   help='The minimum number of ALTs to trigger a detection')
    parser.add_argument('-s','--seed',                 type=int,   help='Use a fixed seed to make the script reproducible')

    args = parser.parse_args()
    
    # Set simulation parameters
    num_pileups     = args.num_pileups 
    depth           = args.mean_total_depth
    error_rate      = args.mean_seq_error_rate
    allele_fraction = args.af
    alt_threshold   = args.alt_threshold
    fv              = args.frac_pileups_with_var

    # Set random seed for reproducibility (optional)
    if args.seed != None:
        np.random.seed(args.seed)

    confusion_matrix, ppa, ppv, specificity = run_simulation_and_gather_metrics(num_pileups, depth, error_rate, allele_fraction, alt_threshold, fv)
    print("Confusion matrix: ")
    print(confusion_matrix)
    print(f"PPA: {ppa}, PPV: {ppv}, Specificity: {specificity}")
    print("="*80)

if __name__ == '__main__':
    main()
