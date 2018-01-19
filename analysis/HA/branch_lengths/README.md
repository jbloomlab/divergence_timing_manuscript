# Branch Length Analysis

This directory contains the analysis of branch lengths for trees optimized under different substitution models.

## Organization

1. [`references`](/.references/): Workflow diagram
2. [`phydms`](./phydms/): output files from `phydms` for the true alignments
3. [`bootstrap`](./bootstrap/): output files from the `phydms` run on the bootstrap alignments.

## Code
The code to optimize the trees and summarize the results is found in the [`snakemake`](http://snakemake.readthedocs.io/en/stable/) file [`Snakefile`](Snakefile).
To run the code, type `snakemake` in *this* directory.
