#!/bin/sh
# properties = {"params": {"r": "3", "prefix": "_WSN_high_0_3", "tree": "../HA/branch_lengths/phydms/WSN_high_0_ExpCM_HA_Doud_prefs_tree.newick"}, "input": [], "output": ["_WSN_high_0_3_simulatedalignment.fasta"], "rule": "simulate_alignment", "local": false, "threads": 1, "cluster": {}, "resources": {}, "jobid": 18}
cd /fh/fast/bloom_j/computational_notebooks/skhilton/2017/divergence_timing_manuscript/analysis/simulations && \
/app/python3/3.4.1/bin/python3.4 -m snakemake _WSN_high_0_3_simulatedalignment.fasta --snakefile /fh/fast/bloom_j/computational_notebooks/skhilton/2017/divergence_timing_manuscript/analysis/simulations/Snakefile \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /fh/fast/bloom_j/computational_notebooks/skhilton/2017/divergence_timing_manuscript/analysis/simulations/.snakemake/tmp.zxlmuba_ --latency-wait 180 \
--benchmark-repeats 1 \
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --quiet --no-hooks --nolock --force-use-threads  --allowed-rules simulate_alignment  && touch "/fh/fast/bloom_j/computational_notebooks/skhilton/2017/divergence_timing_manuscript/analysis/simulations/.snakemake/tmp.zxlmuba_/18.jobfinished" || (touch "/fh/fast/bloom_j/computational_notebooks/skhilton/2017/divergence_timing_manuscript/analysis/simulations/.snakemake/tmp.zxlmuba_/18.jobfailed"; exit 1)

