"""
Python utilities for simulations analysis
Written by Sarah Hilton.
"""
from phydmslib.constants import *
import phydmslib.file_io
from ete3 import Tree
import os
import itertools
from scipy.special import comb
import pandas as pd
import scipy


def createModel(prefs, params):
    """
    This function creates the ExpCM for the `pyvolve` simulation
    input:
        `str`, path to the preferences file
        `str`, path to file with the preference parameters
    """
    prefs = phydmslib.file_io.readPrefs(prefs)
    sites = sorted(prefs.keys())
    prefs = [prefs[r] for r in sites]
    params = pd.read_csv(params, engine="python", sep=" = ", header=None)
    params = dict(zip(params[0], params[1]))
    params["phiT"] = 1 - sum([params[x] for x in params.keys()
                             if x.startswith("phi")])
    phi = scipy.array([params["phi{0}".format(INDEX_TO_NT[x])]
                       for x in range(N_NT)])
    return phydmslib.models.ExpCM(prefs, kappa=params["kappa"],
                                  omega=params["omega"], beta=params["beta"],
                                  mu=0.3, phi=phi, freeparams=['mu'])


def extract_model_name(model_string):
    """
    This function takes a tree output file from phydms and extracts the name of
    the model.
    input:
        `str`, tree file name
    output:
        `str`, model name
        `int`, 1 if "control" model, 0 otherwise
    """
    if "YNGKP" in model_string:
        if "YNGKP_M0" in model_string:
            return "GY94", 0
        elif "YNGKP_M5" in model_string:
            return "GY94_gammaomega", 0
        else:
            raise ValueError("Cannot parse model string {0}"
                             .format(model_string))
    elif "ExpCM" in model_string:
        if "gammaomega" in model_string:
            if "randomized_ExpCM" in model_string:
                return "ExpCM_gammaomega_randomized", 1
            elif "averaged_ExpCM" in model_string:
                return "ExpCM_gammaomega_averaged", 1
            else:
                return "ExpCM_gammaomega", 0
        else:
            if "randomized_ExpCM" in model_string:
                return "ExpCM_randomized", 1
            elif "averaged_ExpCM" in model_string:
                return "ExpCM_averaged", 1
            else:
                return "ExpCM", 0
    elif "RAxML" in model_string:
        return "GTR", 0
    else:
        raise ValueError("Cannot parse model string {0}".format(model_string))


def calc_distances(treePath):
    '''
    This function calculates the distance between all pairs of tips on a given
    tree.
    input:
        treePath: `str`, file path to tree
    output:
        `pandas dataframe`, dataframe summarizing the tree. Columns include the
        names of the two sequences, the identifier for the branch (`seq_id`),
        the two HA groups the sequences come from and the branch length.
    '''
    df = {"sequence1": [], "sequence2": [], "seq_id": [], "distance": []}
    treeName = os.path.basename(treePath)
    with open(treePath) as f:  # workaround for file I/O deprecation in `ete3`
        treeString = f.read()
    t = Tree(treeString)
    leaves = [leaf.name for leaf in t.iter_leaves()]
    for pair in itertools.combinations(leaves, 2):
        seqs = [pair[0], pair[1]]
        seqs.sort()
        df["sequence1"].append(seqs[0])
        df["sequence2"].append(seqs[1])
        df["seq_id"].append("{0}_{1}".format(seqs[0], seqs[1]))
        df["distance"].append(t.get_distance(seqs[0], seqs[1]))
    df = pd.DataFrame(df)
    return df
