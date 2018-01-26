# `phylobayes` analysis

## `phylobayes` mutation-selection model

The code and manual for [`phylobayes`](https://github.com/bayesiancook/pbmpi) can be found [HERE](https://github.com/bayesiancook/pbmpi).

## Code

The analysis is performed by the [`Snakefile`](Snakefile).

There are three steps of the `phylobayes` analysis

1. `pb_mpi`: runs the model
2. `read_mpi`: extracts relevant parameters from the chain
3. `bpcomp`: extracts the consensus tree

After `phylobayes`,  I have other functions which
1. format the preferences
2. calculate the branch lengths
3. normalize the branch lengths

## Data

I ran this protocol on all 7 alignments from HA branch lengths

1. WSN high divergence
2. WSN intermediate divergence
3. WSN low divergence
4. Perth high divergence
5. Perth intermediate divergence
6. Perth low divergence
7. hybrid high divergence 
