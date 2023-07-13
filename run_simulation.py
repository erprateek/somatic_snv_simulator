#!/bin/python
"""
Script Name: run_simulation.py

Monte-Carlo simulation of somatic mutations in NGS data.

This script performs a Monte-Carlo simulation of somatic mutations in 
Next-Generation Sequencing (NGS) data. It takes command-line arguments 
to configure the simulation parameters and generates plots of PPA (Positive 
Predictive Accuracy) and PPV (Positive Predictive Value) as a function of 
the somatic allele fraction.

The script consists of the following functions:
- get_available_memory: Retrieves the available memory on the system and 
returns 60% of it as an integer value.
- calculate_metrics: Calculates various metrics based on the given labels 
and predictions.
- simulate_somatic_mutations: Simulates somatic mutations for a given number 
of pileups and returns the alt reads and mutation presence.
- run_simulation_and_gather_metrics: Runs simulations of somatic mutations 
and gathers metrics based on the results.
- get_chunking_params: Determines the chunking parameters for processing a 
given number of simulations.
- plot_metrics: Plots the PPA and PPV as a function of the somatic 
allele fraction.

Command-line Arguments:
    -k, --num_pileups:          Number of pileups to simulate (default: 10000).
    -D, --mean_total_depth:     Mean of the Poisson distribution, from which depth 
    values are drawn (default: 300).
    -e, --mean_seq_error_rate:  Average sequencing error rate. The number of ALT 
    reads due to noise is drawn from the Binomial distribution (default: 0.005).
    -f, --frac_pileups_with_var: The fraction of pileups with a somatic variant 
    (default: 0.01).
    -a, --af:                   The intrinsic Allele Fraction of somatic variants. 
    The number of ALT reads due to a somatic variant is drawn from the Binomial 
    distribution. Can accept multiple values (default: [0.2]).
    -t, --alt_threshold:        The minimum number of ALTs to trigger a detection 
    (default: 2).
    -s, --seed:                 Use a fixed seed to make the script reproducible.

Note:
    Before running the script, make sure to install the required dependencies listed 
    in the 'requirements.txt' file.

Usage:
    python run_simulation.py [arguments]

Example:
    python run_simulation.py --num_pileups 10000 --mean_total_depth 300 
    --mean_seq_error_rate 0.005 --frac_pileups_with_var 0.01 --af 0.2 
    --alt_threshold 2 --seed 123

Author:
    Prateek Tandon (prateektandon@alumni.cmu.edu)
"""

import argparse
import timeit
import numpy as np
import psutil

def get_available_memory() -> int:
    """
    Retrieves the available memory on the system and returns 60% of it as an 
    integer value.

    Returns:
        int: The available memory as an integer value. Uses only 60% of 
        available memory

    """
    mem = psutil.virtual_memory()
    return int(mem.available*0.6)

def calculate_metrics(labels: np.ndarray, predictions: np.ndarray) -> tuple:
    """
    Calculates various metrics based on the given labels and predictions.

    Args:
        labels (np.ndarray): An array of binary labels.
        predictions (np.ndarray): An array of binary predictions.

    Returns:
        tuple: A tuple containing the true positive, false positive, true 
        negative, and false negative counts.

    """
    labels         = labels.astype('uint64')
    predictions    = predictions.astype('uint64')
    true_positive  = np.sum(labels & predictions)
    false_positive = np.sum(~labels & predictions)
    true_negative  = np.sum(~labels & ~predictions)
    false_negative = np.sum(labels & ~predictions)

    return true_positive, false_positive, true_negative, false_negative

def simulate_somatic_mutations(num_pileups: int,
                               depth: int,
                               error_rate: float,
                               allele_fraction: float,
                               fv: float) -> tuple:
    """
    Simulates somatic mutations for a given number of pileups and returns 
    the alt reads and mutation presence.

    Args:
        num_pileups (int): The total number of pileups to be simulated.
        depth (int): Mean of the poisson distribution from which depth values are drawn.
        error_rate (float): The error rate for simulating sequencing noise.
        allele_fraction (float): The allele fraction for simulating mutation 
        reads.
        fv (float): The fraction of k pileups with a somatic mutation.

    Returns:
        tuple: A tuple containing the alt reads and the presence of mutations.

    """
    # Draw total read depth from a Poisson distribution
    total_depth            = np.random.poisson(depth,
                                               size=num_pileups).astype('uint8')

    # Draw sequencing noise from a Binomial distribution
    noisy_alt_reads        = np.random.binomial(total_depth,
                                                error_rate).astype('uint8')

    # Generate reads having mutation - Using 1 for True and 0 for False
    has_mutation           = np.random.choice([True, False],
                                              size=num_pileups,
                                              p=[fv, 1-fv]).astype('uint8')

    somatic_mutation_reads = np.random.binomial(total_depth,
                                                allele_fraction).astype('uint8')*has_mutation

    total_alt_reads = noisy_alt_reads + somatic_mutation_reads

    return (total_alt_reads, has_mutation)

def run_simulation_and_gather_metrics(num_pileups: int,
                                      depth: int,
                                      error_rate: float,
                                      allele_fraction: float,
                                      alt_threshold: float,
                                      fv: float) -> tuple:
    """
    Runs simulations of somatic mutations and gathers metrics based on the 
    results.

    Args:
        num_pileups (int): The total number of pileups (simulations) to be 
        processed.
        depth (int): The depth parameter for simulating somatic mutations.
        error_rate (float): The error rate for simulating somatic mutations.
        allele_fraction (float): The allele fraction for simulating somatic 
        mutations.
        alt_threshold (float): The alternative allele threshold for determining 
        predictions.
        fv (float): The fraction of k pileups with a somatic variant 

    Returns:
        tuple: A tuple containing the confusion matrix, PPA (Positive Predictive
        Accuracy), 
               PPV (Positive Predictive Value), and specificity.

    """
    num_chunks, chunk_size = get_chunking_params(num_pileups)
    remaining_simulations  = num_pileups % chunk_size
    tp_sum = fp_sum = tn_sum = fn_sum = 0
    for cur_chunk in range(num_chunks):
        print(f"\tProcessing {cur_chunk+1} out of {num_chunks} \
            [{round((cur_chunk+1)*100/(num_chunks), 2)}%]", end='\r')

        # Simulate somatic mutations
        pileups, labels = simulate_somatic_mutations(chunk_size,
                                                    depth,
                                                    error_rate,
                                                    allele_fraction,
                                                    fv)
        predictions = pileups >= alt_threshold
        # tp -> True positives
        # fp -> False positives
        # fn -> False negatives
        # tn -> True negatives
        tp, fp, tn, fn = calculate_metrics(labels, predictions)
        tp_sum += tp
        fp_sum += fp
        tn_sum += tn
        fn_sum += fn
    if remaining_simulations > 0:
        pileups, labels = simulate_somatic_mutations(remaining_simulations,
                                                    depth,
                                                    error_rate,
                                                    allele_fraction,
                                                    fv)
        predictions = pileups >= alt_threshold
        tp, fp, tn, fn = calculate_metrics(labels, predictions)
        tp_sum += tp
        fp_sum += fp
        tn_sum += tn
        fn_sum += fn

    confusion_matrix = np.array([[tn_sum, fp_sum], [fn_sum, tp_sum]])
    ppa              = np.array(tp_sum) / (np.array(tp_sum) + np.array(fn_sum)+1)
    ppv              = np.array(tp_sum) / (np.array(tp_sum) + np.array(fp_sum)+1)
    specificity      = np.array(tn_sum) / (np.array(tn_sum) + np.array(fp_sum)+1)

    # Print the results
    return (confusion_matrix, ppa, ppv, specificity)

def get_chunking_params(total_simulations):
    """
    Determines the chunking parameters for processing a given number of 
    simulations.

    Args:
        total_simulations (int): The total number of simulations to be 
        processed.

    Returns:
        tuple: A tuple containing the number of chunks and the ideal chunk 
        size.

    """
    # Get available memory
    available_memory = get_available_memory()

    # Calculate the ideal chunk size based on available memory
    # Assuming 8 bytes per sample and 3 distributions per simulation
    ideal_chunk_size = min(total_simulations, available_memory // 24)

    # Calculate the number of chunks and remaining simulations
    num_chunks = total_simulations // ideal_chunk_size
    print(f"Available memory: {available_memory}")
    print(f"Automatically assessed ideal chunk size to be: {ideal_chunk_size}")
    print(f"Inferred number of chunks to process: {num_chunks}")
    return (num_chunks, ideal_chunk_size)

def main():
    """
    Monte-Carlo simulation of somatic mutations in NGS data.

    This function performs a Monte-Carlo simulation of somatic mutations in 
    Next-Generation Sequencing (NGS) data. It takes command-line arguments 
    to configure the simulation parameters and prints the resulting confusion 
    matrix, PPA (Positive Predictive Accuracy), PPV (Positive Predictive 
    Value), and specificity.

    Command-line Arguments:
        -k, --num_pileups:           Number of pileups to simulate.
        -D, --mean_total_depth:      Mean of the Poisson distribution, from 
        which depth values are drawn.
        -e, --mean_seq_error_rate:   Average sequencing error rate. The number 
        of ALT reads due to noise is drawn from the Binomial distribution.
        -f, --frac_pileups_with_var: The fraction of pileups with a somatic 
        variant.
        -a, --af:                    The intrinsic Allele Fraction of somatic 
        variants. The number of ALT reads due to a somatic variant is drawn 
        from the Binomial distribution.
        -t, --alt_threshold:         The minimum number of ALTs to trigger a 
        detection.
        -s, --seed:                  Use a fixed seed to make the script 
        reproducible.

    """
    start_time = timeit.default_timer()
    parser = argparse.ArgumentParser("Monte-Carlo simulation of somatic \
                                     mutations in NGS data.")
    parser.add_argument('-k','--num_pileups',          type=int, required=True,
                        help='Number of pileups to simulate')
    parser.add_argument('-D','--mean_total_depth',     type=float, required=True,
                        help="Mean of the Poisson distribution, \
                            from which depth values are drawn")
    parser.add_argument('-e','--mean_seq_error_rate',  type=float, required=True,
                        help="Average sequencing error rate.\
                        The number of ALT reads due to noise is drawn from \
                            BD(D,e)")
    parser.add_argument('-f','--frac_pileups_with_var',type=float, required=True, 
                        help='The fraction of k pileups with a somatic\
                            variant')
    parser.add_argument('-a','--af',                   type=float, required=True,
                        help="""The intrinsic Allele Fraction of somatic\
                            variants.
                        The number of ALT reads due to a somatic variant \
                            is drawn from BD(D,AF)""")
    parser.add_argument('-t','--alt_threshold',        type=int, required=True,
                        help='The minimum number of ALTs to trigger a \
                            detection')
    parser.add_argument('-s','--seed',                 type=int, required=False,
                        help='Use a fixed seed to make the script \
                            reproducible')

    args = parser.parse_args()

    # Set simulation parameters
    num_pileups     = args.num_pileups
    depth           = args.mean_total_depth
    error_rate      = args.mean_seq_error_rate
    allele_fraction = args.af
    alt_threshold   = args.alt_threshold
    fv              = args.frac_pileups_with_var

    # Set random seed for reproducibility (optional)
    if args.seed is not None:
        np.random.seed(args.seed)

    assert num_pileups > 0
    assert depth > 0
    assert error_rate >= 0
    assert allele_fraction >= 0
    assert alt_threshold > 0
    assert fv >= 0

    confusion_matrix, ppa, ppv, specificity = run_simulation_and_gather_metrics(
        num_pileups, 
        depth,
        error_rate, 
        allele_fraction, 
        alt_threshold, 
        fv)
    print("Confusion matrix: ")
    print(confusion_matrix)
    print(f"PPA: {round(ppa*100, 2)}%, \t PPV: {round(ppv*100, 2)}%, \t \
        Specificity: {round(specificity*100, 2)}%")
    print("="*80)
    end_time = timeit.default_timer() - start_time
    print(f"Finished simulating {args.num_pileups} pileups in {end_time} seconds")


if __name__ == '__main__':
    main()
