"""
This script is the scratch script for optimizing a model with different betas.

The script is a modified (simplified) verison of the `phydms` script

SKH 20180122
"""

import sys
import os
import re
import logging
import random
import time
import math
import operator
import multiprocessing
import warnings
warnings.simplefilter('always')
import scipy
import scipy.stats
import statsmodels.sandbox.stats.multicomp
import Bio.Phylo
import phydmslib
from phydmslib.constants import *
import phydmslib.file_io
import phydmslib.parsearguments
import phydmslib.models
import phydmslib.treelikelihood
import pandas as pd

def read_phydms_model_params(param):
    params = pd.read_csv(param, engine="python", sep=" = ", header=None)
    params = dict(zip(params[0], params[1]))
    return params

def create_model_ExpCM(param, prefs, beta):
    """
    This is the basic `ExpCM` from `phydms`

    `omega` is set to 1.0
    """
    params = read_phydms_model_params(param)
    prefs = phydmslib.file_io.readPrefs(prefs)
    sites = sorted(prefs.keys())
    prefs = [prefs[r] for r in sites]
    params["phiT"] = 1 - sum([params[x] for x in params.keys() if
                              x.startswith("phi")])
    phi = scipy.array([params["phi{0}".format(INDEX_TO_NT[x])] for x in
                       range(N_NT)])
    return phydmslib.models.ExpCM(prefs, kappa=params["kappa"],
                                  omega=params["omega"], beta=beta,
                                  mu=0.3, phi=phi, freeparams=['mu'])

def main():
    """Main body of script."""

    alignment = "../../analysis/HA/data/subsample/HA_hybrid_lowH1_0.fasta"
    prefs_fname = "../../analysis/HA/data/references/HA_hybridDoud_prefs.csv"
    params_fname = "../../analysis/HA/branch_lengths/phydms/hybrid_lowH1_0_ExpCM_HA_hybridDoud_prefs_modelparams.txt"
    tree_fname = "../../analysis/HA/branch_lengths/phydms/hybrid_lowH1_0_RAxML_tree.newick"
    # betas = [1.56, 1.34, 1.22]
    betas = [1.00]

    for beta in betas:
        outprefix = "_temp_{0}/".format(beta)
        # create output directory if needed
        outdir = os.path.dirname(outprefix)
        if outdir:
            if not os.path.isdir(outdir):
                if os.path.isfile(outdir):
                    os.remove(outdir)
                os.mkdir(outdir)

        # output files, remove if they already exist
        underscore = '' if outprefix == '/' else '_'
        logfile = '{0}{1}log.log'.format(outprefix, underscore)
        modelparamsfile = '{0}{1}modelparams.txt'.format(outprefix, underscore)
        loglikelihoodfile = '{0}{1}loglikelihood.txt'.format(outprefix, underscore)
        treefile = '{0}{1}tree.newick'.format(outprefix, underscore)
        to_remove = [modelparamsfile, loglikelihoodfile, treefile, logfile,]
        for f in to_remove:
            if os.path.isfile(f):
                os.remove(f)

        # begin execution
        # print some basic information
        random.seed(0)
        # read alignment
        alignment = phydmslib.file_io.ReadCodonAlignment(alignment,
                checknewickvalid=True)
        seqnames = set([head for (head, seq) in alignment])

        # process the substitution model
        model = create_model_ExpCM(params_fname, prefs_fname, beta)

        # read tree
        tree = Bio.Phylo.read(tree_fname, 'newick')
        tipnames = set([clade.name for clade in tree.get_terminals()])
        assert len(tipnames) == tree.count_terminals(), "non-unique tip names?"
        assert tipnames == seqnames, ("Names in alignment do not match those in "
                "tree.\nSequences in alignment but NOT in tree:\n\t{0}\n"
                "Sequences in tree but NOT in alignment:\n\t{1}".format(
                '\n\t'.join(seqnames - tipnames), '\n\t'.join(tipnames - seqnames)))
        tree.root_at_midpoint()
        assert tree.is_bifurcating(), "Tree is not bifurcating: cannot handle"
        nadjusted = 0
        for node in tree.get_terminals() + tree.get_nonterminals():
            if (node.branch_length == None) and (node == tree.root):
                node.branch_length = 1e-06
            elif node.branch_length < 1e-06:
                nadjusted += 1
                node.branch_length = 1e-06

        # set up tree likelihood
        tl = phydmslib.treelikelihood.TreeLikelihood(tree, alignment, model)

        # maximize likelihood
        optimize_brlen = True
        printfunc = None
        maxresult = tl.maximizeLikelihood(optimize_brlen=optimize_brlen,
                approx_grad=False, printfunc=printfunc)
        with open(loglikelihoodfile, 'w') as f:
            f.write('log likelihood = {0:.2f}'.format(tl.loglik))
        params = '\n\t'.join(['{0} = {1:6g}'.format(p, pvalue) for (p, pvalue)
                in sorted(tl.model.paramsReport.items())])
        with open(modelparamsfile, 'w') as f:
            f.write(params.replace('\t', ''))
        Bio.Phylo.write(tl.tree, treefile, 'newick')

if __name__ == '__main__':
    main() # run the script
