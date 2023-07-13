#!/bin/python

import matplotlib.pyplot as plt
import argparse
from run_simulation import run_simulation_and_gather_metrics

def plot_metrics(allele_fractions, ppas, ppvs):
    plt.plot(allele_fractions, ppas, label='PPA')
    plt.plot(allele_fractions, ppvs, label='PPV')
    plt.xlabel('Allele Fraction')
    plt.ylabel('Score')
    plt.legend()
    plt.savefig('Somatic mutation simulation - AF vs PPA and PPV.png', dpi=200)

    
def main():
    parser = argparse.ArgumentParser("Monte-Carlo simulation of somatic mutations in NGS data.")
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

    # Set random seed for reproducibility (optional)
    if args.seed != None:
        np.random.seed(args.seed)

    ppas = []
    ppvs = []

    for af in allele_fractions:
        confusion_matrix, ppa, ppv, specificity = run_simulation_and_gather_metrics(num_pileups, depth, error_rate, af, alt_threshold, fv)
        ppas.append(ppa)
        ppvs.append(ppv)

    plot_metrics(allele_fractions, ppas, ppvs)
    print(list(zip(allele_fractions, ppas, ppvs)))


if __name__ == '__main__':
    main()
