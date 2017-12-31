# Decay to stationary state

The code in this directory computes the "Expected pairwise amino-acid identity at a site" for different models.
This information is used in [Figure 2](../../manuscript/manuscript.pdf).

## Formulas

### Expected pairwise amino-acid identity at a site

The expected pairwise amino-acid identity at a site given a model is
$\sum_a{\sum_{x \in a}{p_{r,x}}\sum_{y \in a}{M\left(t\right)_{xy}}}$

where $a$ are amino acids, $x$ are codons, $p_{r,x}$ is the stationary state of the model, and $M\left(t\right)_{xy}$ is the rate of transition from codon $x$ to codon $y$ after some time $t$.  

### Spielman $\omega_r$

A site-specific variation rate can be inferred from an ExpCM, following [spielman2015the](https://academic.oup.com/mbe/article/32/4/1097/1077799).

This quantity is
$\omega_r = \frac{\sum_x \sum_{y \in N_x}p_{r,x}P_{r,xy}}{\sum_x \sum_{y \in Nx}p_{r,x}Q_{xy}}$

where $x,y$ are codons, $r$ is a site, $p_{r,x}$ is the stationary state, $P_{r,xy}$ is the rate of transition from codon $x$ to codon $y$ at site $r$, $q_{xy}$ is the rate of mutation from codon $x$ to codon $y$ and $N_x$ is the set of codons that are synonymous to $x$ and differ by 1 nucleotide.

## Models

I calculated the expected pairwise amino-acid identity at site for four models

1. `GY94` (`YNGKP_M0`)  
2. `GY94 + $\Gamma\omega$` (`YNGKP_M5`)  
3. `ExpCM`
4. `GY94 + $\omega_r$` (spielman $\omega_r$).

These models are defined by three model parameter files

1. [`YNGKP_M0_modelparams.txt`](YNGKP_M0_modelparams.txt)    
2. [`YNGKP_M5_modelparams.txt`](YNGKP_M5_modelparams.txt)
3. [`ExpCM_HA_Doud_prefs_modelparams.txt`](ExpCM_HA_Doud_prefs_modelparams.txt)

These model parameters come from the `WSN_low_0` run.

The fourth model (`GY94 + $\omega_r$`) is also defined by the parameters in [`YNGKP_M0_modelparams.txt`](YNGKP_M0_modelparams.txt) except for the $\omega$ parameter.

The preferences for the `ExpCM` model are in [`HA_Doud_prefs.csv`](HA_Doud_prefs.csv).

## Script

The code is found in [`decay_to_stationary_plot.py`](decay_to_stationary_plot.py).

## Outputs

1. [`spielman_wr.csv`](spielman_wr.csv): The $\omega_r$ value calculated for each site in the `ExpCM`  

Example:   

site|wr
---|---
1|0.36
2|0.34
3|0.35

2. [`expected_identity_given_time_t.csv`](expected_identity_given_time_t.csv): The expected amino-acid identity at each site given time t. This value is called f in the table.


Example:

Model|Site|Time|f
---|---|---|---|
ExpCM|1|0|1.0
ExpCM|1|1|0.60
ExpCM|1|2|0.41

## Stationary state script

The model feature figure requires  
1. the stationary state of an `ExpCM`
2. the stationary state of a `YNGKP M0`
3. the 4 $\omega_r$ from a `YNGKP M5`
4. the $\omega$ from `YNGKP M0`

The script [`stationary_state.py`](stationary_state.py) extracts the $\omega$ values from the modelparam files and creates a set of "preferences" from the stationary state of an `ExpCM` and `YNGKP M0`. 
