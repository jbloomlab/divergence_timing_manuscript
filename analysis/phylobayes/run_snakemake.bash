#!/bin/bash

snakemake -j 50 --cluster "sbatch -p campus -c 1 -t 10:00:00" --latency-wait 180
