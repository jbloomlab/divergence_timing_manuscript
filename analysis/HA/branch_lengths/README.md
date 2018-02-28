# Branch Length Analysis

This directory contains the analysis of branch lengths for trees optimized under different substitution models.

## Code
The code to optimize the trees and summarize the results is found in the [`snakemake`](http://snakemake.readthedocs.io/en/stable/) file [`Snakefile`](Snakefile).
To run the code, type `snakemake` in *this* directory.

## Organization

[`phydms`](./phydms/): output files from `phydms` for alignments

The [`Snakefile`](Snakefile) creates an `outputs` directory for all of the results. 
