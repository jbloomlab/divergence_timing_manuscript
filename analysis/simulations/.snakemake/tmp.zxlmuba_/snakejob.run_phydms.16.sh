#!/bin/sh
# properties = {"params": {}, "input": ["_WSN_high_0_6_simulatedalignment.fasta"], "output": ["phydms/WSN_high_0_6/modelcomparison.md", "phydms/WSN_high_0_6/"], "rule": "run_phydms", "local": false, "threads": 1, "cluster": {}, "resources": {}, "jobid": 16}
cd /fh/fast/bloom_j/computational_notebooks/skhilton/2017/divergence_timing_manuscript/analysis/simulations && \
/app/python3/3.4.1/bin/python3.4 -m snakemake phydms/WSN_high_0_6/modelcomparison.md phydms/WSN_high_0_6/ --snakefile /fh/fast/bloom_j/computational_notebooks/skhilton/2017/divergence_timing_manuscript/analysis/simulations/Snakefile \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files _WSN_high_0_6_simulatedalignment.fasta /fh/fast/bloom_j/computational_notebooks/skhilton/2017/divergence_timing_manuscript/analysis/simulations/.snakemake/tmp.zxlmuba_ --latency-wait 180 \
--benchmark-repeats 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules run_phydms  && touch "/fh/fast/bloom_j/computational_notebooks/skhilton/2017/divergence_timing_manuscript/analysis/simulations/.snakemake/tmp.zxlmuba_/16.jobfinished" || (touch "/fh/fast/bloom_j/computational_notebooks/skhilton/2017/divergence_timing_manuscript/analysis/simulations/.snakemake/tmp.zxlmuba_/16.jobfailed"; exit 1)

