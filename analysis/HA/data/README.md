# HA data

This directory contains the pre- and post-processed influenza HA sequences.

## Organization
The directory is organized as follows

1. [`./raw/`](./raw/): Raw sequences from the [Influenza Virus Resource Database](https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database).  
2. [`./references/`](./references/): Files need to run the analysis other than the IVR sequences (DMS reference sequences, required sequences, preferences, etc.)
3. [`./aligned/`](./aligned/): Sequences aligned to DMS references in [`./references/`](./references/)    
4. [`./shared/`](.shared/): Aligned sequences which are shared between all three sequences   
5. [`./subsampled/`](./subsampled/): Final alignments which contain a) all target HA groups b) aligned to DMS references c) subsampled to five sequences per HA group with no more than one sequence per year and d) contain all of the required sequences.

## Code

The code to produce the alignments found in [`./aligned/`](./aligned/), [`.shared/`](.shared/), and [`./subsampled/`](./subsampled/) is found in the [`snakemake`](http://snakemake.readthedocs.io/en/stable/) file [`Snakefile`](Snakefile).
To run the code, type `snakemake` in *this* directory.

Here is a workflow of the `Snakefile`
![](references/workflow.png)

### Notes
• I used the standard `phydms_prepalignment` parameters *except* I set a very low `minidentity` threshold of 0.3.  
• For simplicity, I am only going to track [`/.raw/`](./raw/) (as a `.zip`), [`./aligned/`](`./aligned/`), and [`./subsampled/`](./subsampled/). For a clean run, unzip [`/.raw/`](./raw/), remove [`./subsampled/`](./subsampled/), and run the [`Snakefile`](Snakefile). It takes about ~45 minutes for a clean run on the cluster.   
• The alignment step is the longest. To simple generate new replicates or add in new required sequences, unzip [`./aligned`](/aligned/), remove [`./subsampled`], and run [`Snakefile`](Snakefile). This should take only a few minutes.
