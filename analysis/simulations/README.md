# Simulations

## Goal
The purpose of this analysis is to simulate sequences with **site-specific preferences** and re-infer the branch lengths using models from the `ExpCM` and `GY94` families.

## Organization

The code to run the analysis is found in the [`Snakefile`](Snakefile)

## Basic analysis protocol

1. Simulate sequences under the `ExpCM` using [`phydms`](http://jbloomlab.github.io/phydms/) and [`pyvovle`](https://github.com/sjspielman/pyvolve).
2. Re-infer the branches using [`phydms_comprehensive`](http://jbloomlab.github.io/phydms/).
3. Calculate the branches from every pair of sequences.
4. Calculate the mean branch length for a given model across all of the simulations.
