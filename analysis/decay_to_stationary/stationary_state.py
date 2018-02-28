"""
The purpose of this script is to create "preferences" for the stationary state of the
YNGKP M0, YNGKP M5, and the ExpCM model.

SKH 20171215
"""

import pandas as pd
from phydmslib.constants import *
import phydmslib.models
import phydmslib.file_io
import subprocess
import os

def read_phydms_model_params(param):
    params = pd.read_csv(param, engine="python", sep= " = ", header=None)
    params = dict(zip(params[0], params[1]))
    return params

def create_model_YNGKP_M0(param, nsites):
    """
    This is the basic `YNGKP_M0` from `phydms`
    """
    params = read_phydms_model_params(param)
    e_pw = [[] for x in range(3)]
    for i in range(3):
        params["phi{0}T".format(i)] = 1 - sum([params[x] for x in params.keys() if x.startswith("phi{0}".format(i))])
        e_pw[i] = [params["phi{0}{1}".format(i, INDEX_TO_NT[x])] for x in range(N_NT)]
    e_pw = scipy.array(e_pw)
    return phydmslib.models.YNGKP_M0(e_pw, nsites, kappa=params["kappa"], omega=params["omega"], mu=0.3,
            freeparams=['mu'])

def create_model_ExpCM(param, prefs):
    """
    This is the basic `ExpCM` from `phydms`

    `omega` is set to 1.0
    """
    params = read_phydms_model_params(param)
    prefs = phydmslib.file_io.readPrefs(prefs)
    sites = sorted(prefs.keys())
    prefs = [prefs[r] for r in sites]
    params["phiT"] = 1 - sum([params[x] for x in params.keys() if x.startswith("phi")])
    phi = scipy.array([params["phi{0}".format(INDEX_TO_NT[x])] for x in range(N_NT)])
    return phydmslib.models.ExpCM(prefs, kappa=params["kappa"], omega=params["omega"],
        beta=params["beta"], mu=0.3, phi=phi, freeparams=['mu'])

def get_stationarystate(model,start,length):
    amino_acids = list(AA_TO_INDEX.keys())
    prefs = {}
    for amino_acid in amino_acids:
        prefs[amino_acid] = []
    prefs["site"] = []
    amino_acids.sort()
    final_cols = ["site"] + amino_acids

    # extract the stationary state
    for r in range(start, start + length):
        prefs["site"].append(r)
        for amino_acid in amino_acids:
            amino_acid_index = AA_TO_INDEX[amino_acid]
            codon_indices = scipy.where(CODON_TO_AA==amino_acid_index)[0]
            prefs[amino_acid].append(model.stationarystate[r][codon_indices].sum())
    prefs = pd.DataFrame(prefs)
    prefs = prefs[final_cols]
    prefs["site"] = [x+1 for x in range(len(prefs))]
    return prefs

def main():
    # set up params for runs
    start_site = 117
    number_of_sites = 15
    phydms_dir = "../HA/branch_lengths/phydms/"
    prefs_dir = "../HA/data/references/"
    YNGKP_M0_modelparams_fname = "{0}hybrid_lowH1_0_YNGKP_M0_modelparams.txt".format(phydms_dir)
    YNGKP_M5_modelaprams_fname = "{0}hybrid_lowH1_0_YNGKP_M5_modelparams.txt".format(phydms_dir)
    ExpCM_modelparams_fname = "{0}hybrid_lowH1_0_ExpCM_HA_hybridDoud_prefs_modelparams.txt".format(phydms_dir)
    Doud_prefs_fname = "{0}HA_hybridDoud_prefs.csv".format(prefs_dir)
    if not os.path.isdir("outputs"):
        os.makedirs("outputs")

    # ExpCM
    model = create_model_ExpCM(ExpCM_modelparams_fname, Doud_prefs_fname)
    prefs = get_stationarystate(model, start_site, number_of_sites)
    prefs.to_csv("outputs/ExpCM_stationarystate.csv", index=False)

    # YNGKP M0
    model = create_model_YNGKP_M0(YNGKP_M0_modelparams_fname, 565)
    prefs = get_stationarystate(model, start_site, number_of_sites)
    prefs.to_csv("outputs/YNGKP_M0_stationarystate.csv", index=False)

    # extract the wr values from YNGKP M5 and YNGKP
    df = pd.read_csv(YNGKP_M5_modelaprams_fname, sep=" = ", engine="python", header=None)
    alpha_omega = df[df[0] == "alpha_omega"][1].iloc[0]
    beta_omega = df[df[0] == "beta_omega"][1].iloc[0]
    wr = list(phydmslib.models.DiscreteGamma(alpha_omega, beta_omega, 4))
    wr.extend([model.omega, alpha_omega, beta_omega])
    value_type = ["omega","omega", "omega", "omega", "omega", "gamma", "gamma"]
    model_list = ["M5", "M5", "M5", "M5", "M0", "alpha_omega", "beta_omega"]
    df = pd.DataFrame({"model":model_list, "omega":wr, "type":value_type})
    df.to_csv("outputs/wr.csv", index=False)

if __name__ == '__main__':
    main()
