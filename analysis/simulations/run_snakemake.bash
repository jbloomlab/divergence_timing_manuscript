#!/bin/bash

snakemake -j 50 --cluster "sbatch -p largenode --mem=300000 -c 6" --latency-wait 180
