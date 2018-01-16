#!/bin/sh
# properties = {"local": false, "rule": "run_phydms", "resources": {}, "threads": 1, "cluster": {}, "input": ["_WSN_high_0_6_modelcomparison.md_simulatedalignment.fasta"], "jobid": 20, "params": {}, "output": ["phydms/WSN_high_0_6_modelcomparison.md/modelcomparison.md", "phydms/WSN_high_0_6_modelcomparison.md"]}
cd /fh/fast/bloom_j/computational_notebooks/skhilton/2017/divergence_timing_manuscript/analysis/simulations && \
/app/python3/3.4.1/bin/python3.4 -m snakemake phydms/WSN_high_0_6_modelcomparison.md/modelcomparison.md phydms/WSN_high_0_6_modelcomparison.md --snakefile /fh/fast/bloom_j/computational_notebooks/skhilton/2017/divergence_timing_manuscript/analysis/simulations/Snakefile \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files _WSN_high_0_6_modelcomparison.md_simulatedalignment.fasta /fh/fast/bloom_j/computational_notebooks/skhilton/2017/divergence_timing_manuscript/analysis/simulations/.snakemake/tmp.ih1tumzp --latency-wait 180 \
--benchmark-repeats 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules run_phydms  && touch "/fh/fast/bloom_j/computational_notebooks/skhilton/2017/divergence_timing_manuscript/analysis/simulations/.snakemake/tmp.ih1tumzp/20.jobfinished" || (touch "/fh/fast/bloom_j/computational_notebooks/skhilton/2017/divergence_timing_manuscript/analysis/simulations/.snakemake/tmp.ih1tumzp/20.jobfailed"; exit 1)

