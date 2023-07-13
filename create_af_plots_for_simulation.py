#!/bin/python
"""
Script Name: create_af_plots_for_simulation.py

Description:
    Monte-Carlo simulation of somatic mutations in NGS data - AF plots version

Command-line Arguments:
    -k, --num_pileups:           Number of pileups to simulate 
    (default: 10000).
    -D, --mean_total_depth:      Mean of the Poisson distribution, 
    from which depth values are drawn (default: 300).
    -e, --mean_seq_error_rate:   Average sequencing error rate. The 
    number of ALT reads due to noise is drawn from the Binomial 
    distribution (default: 0.005).
    -f, --frac_pileups_with_var: The fraction of pileups with a somatic 
    variant (default: 0.01).
    -a, --af:                    The intrinsic Allele Fraction of somatic 
    variants. The number of ALT reads due to a somatic variant is drawn 
    from the Binomial distribution. Can accept multiple values (default: 
    [0.2]).
    -t, --alt_threshold:         The minimum number of ALTs to trigger a 
    detection (default: 2).
    -s, --seed:                  Use a fixed seed to make the script 
    reproducible.
    
Author:
    Prateek Tandon (prateektandon@alumni.cmu.edu)

Dependencies:
    - matplotlib
    - numpy
    - run_simulation

Usage:
    python create_af_plots_for_simulation.py [arguments]

"""
import argparse
import matplotlib.pyplot as plt
import numpy as np
from run_simulation import run_simulation_and_gather_metrics

def plot_metrics(allele_fractions: np.ndarray, 
                 ppas: np.ndarray, ppvs: np.ndarray) -> None:
    """
    Plots the PPA (Positive Predictive Accuracy) and PPV (Positive 
    Predictive Value) as a function of the somatic allele fraction.

    Args:
        allele_fractions (np.ndarray): An array of somatic allele fractions.
        ppas (np.ndarray): An array of corresponding PPA values.
        ppvs (np.ndarray): An array of corresponding PPV values.

    Returns:
        None

    """
    plt.plot(allele_fractions, ppas, label='PPA', marker='o')
    plt.plot(allele_fractions, ppvs, label='PPV', marker='x')
    plt.xlabel('Allele Fraction')
    plt.ylabel('Score')
    plt.title('PPA and PPV as a function of the somatic AF')
    plt.legend()
    plt.savefig('Somatic mutation simulation - AF vs PPA and PPV.png', dpi=200)
   
def main():
    """
    Monte-Carlo simulation of somatic mutations in NGS data - AF plots version.

    This function performs a Monte-Carlo simulation of somatic mutations in 
    Next-Generation Sequencing (NGS) data. It takes command-line arguments to 
    configure the simulation parameters and generates plots of PPA (Positive 
    Predictive Accuracy) and PPV (Positive Predictive Value) as a function 
    of the somatic allele fraction.

    Command-line Arguments:
        -k, --num_pileups:          Number of pileups to simulate (default: 10000).
        -D, --mean_total_depth:     Mean of the Poisson distribution, from which depth values are drawn (default: 300).
        -e, --mean_seq_error_rate:  Average sequencing error rate. The number of ALT reads due to noise is drawn from the Binomial distribution (default: 0.005).
        -f, --frac_pileups_with_var: The fraction of pileups with a somatic variant (default: 0.01).
        -a, --af:                   The intrinsic Allele Fraction of somatic variants. The number of ALT reads due to a somatic variant is drawn from the Binomial distribution. Can accept multiple values (default: [0.2]).
        -t, --alt_threshold:        The minimum number of ALTs to trigger a detection (default: 2).
        -s, --seed:                 Use a fixed seed to make the script reproducible.

    """
    
    parser = argparse.ArgumentParser("Monte-Carlo simulation of somatic mutations in NGS data - AF plots version")
    parser.add_argument('-k','--num_pileups',          type=int,   default=10000, help='Number of pileups to simulate')
    parser.add_argument('-D','--mean_total_depth',     type=float, default=300,   help='Mean of the Poisson distribution, from which depth values are drawn')
    parser.add_argument('-e','--mean_seq_error_rate',  type=float, default=0.005, help='Average sequencing error rate. The number of ALT reads due to noise is drawn from BD(D,e)')
    parser.add_argument('-f','--frac_pileups_with_var',type=float, default=0.01,  help='The fraction of k pileups with a somatic variant')
    parser.add_argument('-a','--af',                   type=float, default=[0.2], help='The intrinsic Allele Fraction of somatic variants. The number of ALT reads due to a somatic variant is drawn from BD(D,AF)', nargs='+')
    parser.add_argument('-t','--alt_threshold',        type=int,   default=2,     help='The minimum number of ALTs to trigger a detection')
    parser.add_argument('-s','--seed',                 type=int,   help='Use a fixed seed to make the script reproducible')

    args = parser.parse_args()

    # Set simulation parameters
    num_pileups      = args.num_pileups
    depth            = args.mean_total_depth
    error_rate       = args.mean_seq_error_rate
    allele_fractions = args.af
    alt_threshold    = args.alt_threshold
    fv               = args.frac_pileups_with_var

    # Set random seed for reproducibility
    if args.seed is not None:
        np.random.seed(args.seed)

    ppas = []
    ppvs = []

    for af in allele_fractions:
        _, ppa, ppv, _ = \
        run_simulation_and_gather_metrics(num_pileups, 
                                          depth, 
                                          error_rate, 
                                          af, 
                                          alt_threshold, 
                                          fv)
        ppas.append(ppa)
        ppvs.append(ppv)

    plot_metrics(allele_fractions, ppas, ppvs)
    print(list(zip(allele_fractions, ppas, ppvs)))


if __name__ == '__main__':
    main()
