"""
The purpose of this script is to calculate how long it takes a model to reach
an expected pairwise identity of 50%.

SKH 20171215
"""

import pandas as pd
from phydmslib.constants import *
import phydmslib.models
import phydmslib.file_io


## Define the models
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

def create_model_YNGKP_M5(param, nsites):
    """
    This is the basic `YNGKP_M5` from `phydms` with 4 categories.
    """
    params = read_phydms_model_params(param)
    e_pw = [[] for x in range(3)]
    for i in range(3):
        params["phi{0}T".format(i)] = 1 - sum([params[x] for x in params.keys() if x.startswith("phi{0}".format(i))])
        e_pw[i] = [params["phi{0}{1}".format(i, INDEX_TO_NT[x])] for x in range(N_NT)]
    e_pw = scipy.array(e_pw)
    m =  phydmslib.models.YNGKP_M0(e_pw, nsites, kappa=params["kappa"], mu=0.3,
            freeparams=['mu', 'omega'])
    return phydmslib.models.GammaDistributedOmegaModel(m, ncats=4, alpha_lambda=params["alpha_omega"], beta_lambda=params["beta_omega"])

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

def create_model_YNGKP_wr(param,spielman_wr):
    """
    This function creates a `YNGKP_M0` with `omega` equal to the value of
    `spielman_wr`.
    """
    params = read_phydms_model_params(param)
    e_pw = [[] for x in range(3)]
    for i in range(3):
        params["phi{0}T".format(i)] = 1 - sum([params[x] for x in params.keys() if x.startswith("phi{0}".format(i))])
        e_pw[i] = [params["phi{0}{1}".format(i, INDEX_TO_NT[x])] for x in range(N_NT)]
    e_pw = scipy.array(e_pw)
    return phydmslib.models.YNGKP_M0(e_pw, 1, kappa=params["kappa"], omega=spielman_wr, mu=0.3,
            freeparams=['mu'])

## Calculations
def f_calculation(x, pr, Mt, model_name):
    """
    This function calcuates the value of `f` (the y-axis) for some amino-acid `x`,
    site `r`, and time `t`.

    The inputs for this function are the stationary state subsetted to site `r`
    and the transition matrix (M) subsetted to site `r` and time `t`.

    If the model is `YNGKP_M5`, then the input is a list of stationary states
    and transition matrices the length of the number of categories. The `f` value
    will be the average `f` value for the categories.

    >>> print(f_calculation(8, scipy.array([0.1,0.2,0.3,0.4]), scipy.array([[1,-100,3,-100],[-100,-100,-100,-100],[9,-100,11,-100],[-100,-100,-100,-100]]), "ExpCM"))
    6.4
    """
    target_codons = scipy.where(CODON_TO_AA==x)[0] # codons which encode amino-acid `x`
    if model_name == "GY94 + Gr": # need to take the average `f` over categories `k`
        prx = [pr[k][target_codons] for k in range(len(pr))] # stationary state
        Mtxy = [Mt[k][target_codons][:,target_codons] for k in range(len(Mt))] # transition to syn. codons
        Mtxy = [scipy.sum(Mtxy[k], axis=1) for k in range(len(Mt))] # sum the transition to syn codons
        # multiply sum of transition to syn codons by stationary state and sum for each category k
        f_rtx = [(Mtxy[k] * prx[k]).sum() for k in range(len(Mt))]
        return sum(f_rtx)/len(Mt) # return the average `f` for the categories `k`
    elif model_name in ["ExpCM", "GY94", "GY94 + wr"]:
        prx = pr[target_codons] # stationary state
        Mtxy = Mt[target_codons][:,target_codons] # transition to syn codons
        Mtxy = scipy.sum(Mtxy, axis=1) # sum of transition to syn codons
        # multiply sum of the transition to syn codons by the stationary state and sum
        f_rtx = (prx * Mtxy).sum()
        # f_rtx = (pr[target_codons] * scipy.sum(Mt[target_codons][:,target_codons], axis=1)).sum()
    else:
        raise ValueError("Cannot handle model {0}".format(model_name)) # wrong model
    return f_rtx

def get_pr(r,model,model_name):
    """
    The purpose of this function is to return the stationary state subsetted to
    site `r` correctly for each model.

    It returns a (61,61) array for every model except `YNGKP_M5` which it returns
    [(61,61) for k in n.cats].
    """
    if model_name == "GY94 + Gr":
        return [model.stationarystate(k)[r] for k in range(model.ncats)]
    elif model_name in ["ExpCM", "GY94", "GY94 + wr"]:
        return model.stationarystate[r]
    else:
        raise ValueError("Cannot handel model {0}".format(model_name))

def get_Mrt(r,t,model,model_name):
    """
    The purpose of this function is to return the transition matrix subsetted to
    site `r` and time `t`correctly for each model.

    This function also uses model.branchScale to correct `t` for each model.

    It returns a (61,61) array for every model except `YNGKP_M5` which it returns
    [(61,61) for k in n.cats].
    """
    if model_name == "GY94 + Gr":
        return [model.M(k,float(t/model.branchScale))[r] for k in range(model.ncats)]
    elif model_name in ["ExpCM", "GY94", "GY94 + wr"]:
        return model.M(float(t/model.branchScale))[r]
    else:
        raise ValueError("Cannot handel model {0}".format(model_name))

def main():
    target_sites = [90, 51, 23, 490, 81]
    target_divs = [0.95, 0.75, 0.5]
    df = pd.read_csv("expected_identity_given_time_t_long.csv")
    print(len(df))
    df = df[df["Site"].isin(target_sites)]
    print(len(df))

    final = {"Model":[], "Site":[], "Time":[], "f":[], "target_div":[]}
    for target_div in target_divs:
        for name, group in df.groupby(["Model","Site"]):
            time = group[group["f"] >= target_div]["Time"].max()
            final["Model"].append(name[0])
            final["Site"].append(name[1])
            final["Time"].append(time)
            final["f"].append(group[group["f"] >= target_div]["f"].min())
            final["target_div"].append(target_div)
    final = pd.DataFrame(final)
    print(final)

if __name__ == '__main__':
    main()
