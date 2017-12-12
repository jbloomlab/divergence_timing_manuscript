# Reference Files

This directory contains the reference files need to run the alignment creation pipeline.

## Deep Mutational Scanning References
This directory contains the data from the deep mutational scans of HA.
This includes the scan performed by Mike Doud and Jesse Bloom ([doud2016accurate](http://www.mdpi.com/1999-4915/8/6/155)) in the HA H1 strain WSN and the scan performed by Juhye Lee and Jesse Bloom (lee2017Perth) in the HA H3 strain Perth/09. This directory includes **the wildtype sequence for each scan** as well as the **site-specific amino-acid preferences measured by each scan**.

### Wildtype sequences
[`WSN_HA_reference.fasta`](WSN_HA_reference.fasta) and `Perth_HA_reference.fasta`(Perth_HA_reference.fasta) are the nucleotide sequences used by [doud2016accurate](http://www.mdpi.com/1999-4915/8/6/155) and [lee2016Perth]() to perform the DMS.

[`WSN_hybrid_reference.fasta`](WSN_hybrid_reference.fasta) is the WSN reference sequence but only contains the sites *shared* between WSN and Perth.


### Amino-acid preferences

#### [`HA_Lee_prefs.csv`](HA_Lee_prefs.csv)

These measurements are from Juhye's `dms_tools` (not `dms_tools2`) analysis.

#### [`HA_Doud_prefs`](HA_Doud_prefs.csv)

These measurements are from  `/fh/fast/bloom_j/computational_notebooks/mdoud/2016/2016_WSN_HA_analysis/preferences/2016-WSNHA_merged_prefs_renumbered.txt`
Unlike the preference file in Mike's paper, this file is *just* Mike's 2016 data and not the average of Doud2016 and Bargavi's data.

#### Hybrid preferences

The mapping of shared sites to sequential sites can be found in [`hybrid_numbering.csv`](hybrid_numbering.csv).

To make [`HA_average_prefs.csv`](HA_average_prefs.csv) I calculated $\pi_{r,A\left(x\right)} = \frac{\pi_{Doud,r,A\left(x\right)} + \pi_{Lee,r,A\left(x\right)}}{2}$ for $r \in \{\rm{shared \ sites}\}$.

## Required Sequences

The [`required_seqs.csv`](required_seqs.csv) file contains contains sequences which I want to be in *every* final alignment in [`subsample`](../subsample/). This includes the DMS reference strains and "landmark" sequences such as the H3 equine clade, early and late H1 sequences, etc.

## Workflow

The workflow diagram ([`workflow.png`](workflow.png)) and the `powerpoint`([`workflow.pptx`](workflow.pptx)) to create the diagram.  
