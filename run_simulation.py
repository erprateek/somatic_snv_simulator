import numpy as np

def simulate_somatic_mutations(num_pileups, depth, error_rate, allele_fraction, alt_threshold, fv):
    # Simulate total read depth from a Poisson distribution
    total_depth = np.random.poisson(depth, size=num_pileups)

    # Simulate sequencing noise from a Binomial distribution
    noise_reads = np.random.binomial(total_depth, error_rate)

    # Simulate somatic mutation (if any)
    has_mutation = np.random.choice([True, False], size=num_pileups, p=[allele_fraction, 1-allele_fraction])
    #mutation_reads = np.random.binomial(total_depth, allele_fraction) * has_mutation
    mutation_reads = np.random.binomial(total_depth, allele_fraction) * (np.random.rand(num_pileups) < fv) * has_mutation


    # Calculate ALT reads (sequencing noise + somatic mutation)
    alt_reads = noise_reads + mutation_reads

    # Determine if somatic variant is detected
    is_positive = alt_reads >= alt_threshold

    return alt_reads, is_positive

# Set simulation parameters
num_pileups = 10000  # 
depth = 100
error_rate = 0.005
allele_fraction = 0.05
alt_threshold = 2
fv = 0.01

# Set random seed for reproducibility (optional)
np.random.seed(42)

# Simulate somatic mutations
pileups, labels = simulate_somatic_mutations(num_pileups, 
depth, 
error_rate, 
allele_fraction, 
alt_threshold, 
fv)

# Print the results
for i in range(num_pileups):
    print("Pileup:", pileups[i], "Label:", labels[i])
