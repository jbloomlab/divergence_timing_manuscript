#!/bin/bash

snakemake -j 50 --cluster "sbatch -p campus -c 1" --latency-wait 180
