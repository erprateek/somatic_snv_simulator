import numpy as np
import argparse
import matplotlib.pyplot as plt
import psutil

def get_available_memory():
    mem = psutil.virtual_memory()
    return int(mem.available*0.8)

def simulate_somatic_mutations_chunk(chunk_size, depth, error_rate, allele_fraction, alt_threshold, fv):
    # Simulate total read depth from a Poisson distribution
    total_depth = np.random.poisson(depth, size=chunk_size)

    # Simulate sequencing noise from a Binomial distribution
    noise_reads = np.random.binomial(total_depth, error_rate)

    # Simulate somatic mutation (if any)
    has_mutation = np.random.choice([True, False], size=chunk_size, p=[allele_fraction, 1-allele_fraction])
    mutation_reads = np.random.binomial(total_depth, allele_fraction) * (np.random.rand(chunk_size) < fv) * has_mutation

    # Calculate ALT reads (sequencing noise + somatic mutation)
    alt_reads = noise_reads + mutation_reads

    # Determine if somatic variant is detected
    is_positive = alt_reads >= alt_threshold

    return alt_reads, is_positive

def calculate_metrics(labels, predictions):
    true_positive = np.sum(labels & predictions)
    false_positive = np.sum(~labels & predictions)
    true_negative = np.sum(~labels & ~predictions)
    false_negative = np.sum(labels & ~predictions)

    return true_positive, false_positive, true_negative, false_negative

def plot_metrics(allele_fractions, ppas, ppvs):
    plt.plot(allele_fractions, ppas, label='PPA')
    plt.plot(allele_fractions, ppvs, label='PPV')
    plt.xlabel('Allele Fraction')
    plt.ylabel('Score')
    plt.legend()
    plt.show()

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Monte-Carlo simulation of somatic mutations in NGS data')
parser.add_argument('--num_simulations', type=int, default=1000, help='Number of simulations')
parser.add_argument('--depth', type=int, default=50, help='Total read depth')
parser.add_argument('--error_rate', type=float, default=0.05, help='Sequencing error rate')
parser.add_argument('--alt_threshold', type=int, default=3, help='ALT threshold')
parser.add_argument('--allele_fractions', nargs='+', type=float, default=[0.1], help='Allele fractions to test')
parser.add_argument('--fv', type=float, default=0.5, help='Fraction of k pileups with a somatic variant')

args = parser.parse_args()

# Set random seed for reproducibility (optional)
np.random.seed(42)

# Get available memory
available_memory = get_available_memory()

# Calculate the ideal chunk size based on available memory
total_simulations = args.num_simulations
ideal_chunk_size = min(total_simulations, available_memory // 24)  # Assuming 8 bytes per sample

# Calculate the number of chunks and remaining simulations
num_chunks = total_simulations // ideal_chunk_size
remaining_simulations = total_simulations % ideal_chunk_size
print(f"Available memory: {available_memory}")
print(f"Automatically assessed ideal chunk size to be: {ideal_chunk_size}")
print(f"Inferred number of chunks to process: {num_chunks}")

# Initialize variables
true_positive_total = []
false_positive_total = []
true_negative_total = []
false_negative_total = []

# Perform simulations for each allele fraction
for allele_fraction in args.allele_fractions:
    print(f"Processing AF: {allele_fraction}")
    # Initialize metrics for the current allele fraction
    true_positive_sum = 0
    false_positive_sum = 0
    true_negative_sum = 0
    false_negative_sum = 0

    # Perform simulations
    for i in range(num_chunks):
        print(f"\tProcessing {i+1}th chunk out of {num_chunks}")
        # Simulate somatic mutations for the chunk
        pileups, labels      = simulate_somatic_mutations_chunk(ideal_chunk_size,
                                                                args.depth,
                                                                args.error_rate,
                                                                allele_fraction,
                                                                args.alt_threshold,
                                                                args.fv)

        # Perform predictions (using alt_threshold as the decision threshold)
        predictions = pileups >= args.alt_threshold

        # Calculate metrics for the current simulation
        true_positive, false_positive, true_negative, false_negative = calculate_metrics(labels, predictions)

        # Accumulate metrics
        true_positive_sum   += true_positive
        false_positive_sum  += false_positive
        true_negative_sum   += true_negative
        false_negative_sum  += false_negative

    # Handle remaining simulations (if any)
    if remaining_simulations > 0:
        pileups, labels = simulate_somatic_mutations_chunk(remaining_simulations, depth, error_rate, allele_fraction, alt_threshold)
        predictions = pileups >= alt_threshold
        true_positive, false_positive, true_negative, false_negative = calculate_metrics(labels, predictions)
        
        true_positive_total += true_positive
        false_positive_total += false_positive
        true_negative_total += true_negative
        false_negative_total += false_negative


    # Accumulate metrics for the current allele fraction
    true_positive_total.append(true_positive_sum)
    false_positive_total.append(false_positive_sum)
    true_negative_total.append(true_negative_sum)
    false_negative_total.append(false_negative_sum)

# Calculate overall metrics
#ppas = np.array(true_positive_total) / (np.array(true_positive_total) + np.array(false_negative_total))
#ppvs = np.array(true_positive_total) / (np.array(true_positive_total) + np.array(false_positive_total))

total_samples = total_simulations
confusion_matrix = np.array([[true_negative_total, false_positive_total], [false_negative_total, true_positive_total]])
ppas = np.array(true_positive_total) / (np.array(true_positive_total) + np.array(false_negative_total))
ppvs = np.array(true_positive_total) / (np.array(true_positive_total) + np.array(false_positive_total))
specificities = np.array(true_negative_total) / (np.array(true_negative_total) + np.array(false_positive_total))

print("Confusion Matrix:")
print(confusion_matrix)
print("PPA:", ppas)
print("PPV:", ppvs)
print("Specificity:", specificities)

# Plot metrics
plot_metrics(args.allele_fractions, ppas, ppvs)
