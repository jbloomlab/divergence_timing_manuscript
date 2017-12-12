#!/bin/bash

snakemake -j 50 --cluster "sbatch -p campus" --latency-wait 180
