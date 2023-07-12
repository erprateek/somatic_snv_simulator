# Somatic SNV Simulation

Write code to implement a Monte-Carlo simulation of somatic mutations in NGS data.


We will simulate a large number of read pileups, as a list of basecalls that are either the REF allele or an ALT allele.  The Total Read Depth for each pileup is a random integer, drawn from a Poisson distribution.  For the purposes of this simulation, we do not care what the nucleotide alleles are, we will only record REF or ALT for each basecall (which could be encoded with 0 and 1, or False and True, for example).


There is the possibility of sequencing noise in every pileup. Simulate this with a random draw from a Binomial Distribution BD(D,e) where D is the total Depth, and e is the average sequencing error rate.  This will tell us how many of the reads in the pileup need to be recorded as ALT, due to sequencing error.


There is also the possibility that some pileups will contain a somatic mutation. For a random fraction of pileups, we will simulate a somatic mutation, as a random draw from a Binomial Distribution BD(D,AF) where D is the total Depth, and AF is the Allele Fraction of the somatic variant. This will tell us how many reads in the pileup need to be recorded as ALT, due to the presence of a somatic variant.  Whenever a pileup is selected to have a somatic variant, we label that pileup as a Positive for having a somatic variant (even if the binomial draw adds Zero ALT reads!). Likewise, if a pileup is NOT selected to have a somatic variant, we label it as a Negative.


For each simulated pileup, we will determine whether a somatic variant is detected by simply applying a minimum threshold on the number of ALT reads in the pileup.


Parameter table:
----------------
|Parameter|Description|Note|
|---------|-----------|----|
|k|Number of pileups to simulate||
|D|Average Total Depth in a pileup|The depth of any one pileup is drawn from a Poisson Distribution with mean D|
|e|Average sequencing error rate|The number of ALT reads due to noise is drawn from BD(D,e)|
|Fv|The fraction of k pileups with a somatic variant||
|AF|The intrinsic Allele Fraction of somatic variants|The number of ALT reads due to a somatic variant is drawn from BD(D,AF)|
|DT|The minimum number of ALTs to trigger a detection||

Implement the framework for this simulation in Python or R.


1) Set parameters k=10000, D=100, e=0.005, Fv=0.01, AF=0.05, and DT=2. 
Show the confusion matrix, and report PPA, PPV, and Specificity


2) Set parameters k=10000, D=300, e=0.005, Fv=0.01, AF=0.05, and DT=2. 
Show the confusion matrix, and report PPA, PPV, and Specificity


3) Run six simulations where all parameters are held fixed except AF, which is set to:
(0.2, 0.1, 0.05, 0.02, 0.01, 0.005). 
What do PPA and PPV look like as a function of the somatic AF?


4) BONUS
Can you make the framework efficient enough to run a simulation with k=3-billion in a reasonable span of wallclock time?
