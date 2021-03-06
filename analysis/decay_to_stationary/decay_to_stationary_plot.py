"""
The purpose of this script is to perform the calculations to make the
"decay to stationary state" plot.

This includes calculating the "spielman w_r" rates of non-syn change from the
`ExpCM` applying these site-specific values to the `YNGKP_M0`.

To run doctest, type:
`python3 -m doctest -v test_doctest.py`

SKH 20170910
"""

import pandas as pd
from phydmslib.constants import *
import phydmslib.models
import phydmslib.file_io
import os


# Define the models
def read_phydms_model_params(param):
    params = pd.read_csv(param, engine="python", sep=" = ", header=None)
    params = dict(zip(params[0], params[1]))
    return params


def create_model_YNGKP_M0(param, nsites):
    """
    This is the basic `YNGKP_M0` from `phydms`
    """
    params = read_phydms_model_params(param)
    e_pw = [[] for x in range(3)]
    for i in range(3):
        params["phi{0}T".format(i)] = 1 - sum([params[x] for x in params.keys()
                                               if x.startswith("phi{0}"
                                               .format(i))])
        e_pw[i] = [params["phi{0}{1}".format(i, INDEX_TO_NT[x])] for x in
                   range(N_NT)]
    e_pw = scipy.array(e_pw)
    return phydmslib.models.YNGKP_M0(e_pw, nsites, kappa=params["kappa"],
                                     omega=params["omega"], mu=0.3,
                                     freeparams=['mu'])


def create_model_YNGKP_M5(param, nsites):
    """
    This is the basic `YNGKP_M5` from `phydms` with 4 categories.
    """
    params = read_phydms_model_params(param)
    e_pw = [[] for x in range(3)]
    for i in range(3):
        params["phi{0}T".format(i)] = 1 - sum([params[x] for x in params.keys()
                                               if x.startswith("phi{0}"
                                               .format(i))])
        e_pw[i] = [params["phi{0}{1}".format(i, INDEX_TO_NT[x])] for x in
                   range(N_NT)]
    e_pw = scipy.array(e_pw)
    m = phydmslib.models.YNGKP_M0(e_pw, nsites, kappa=params["kappa"], mu=0.3,
                                  freeparams=['mu', 'omega'])
    return phydmslib.models.GammaDistributedOmegaModel(m, ncats=4,
                                                       alpha_lambda=params["alpha_omega"],
                                                       beta_lambda=params["beta_omega"])


def create_model_ExpCM(param, prefs):
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
                                  omega=params["omega"], beta=params["beta"],
                                  mu=0.3, phi=phi, freeparams=['mu'])


def create_model_YNGKP_wr(param, spielman_wr):
    """
    This function creates a `YNGKP_M0` with `omega` equal to the value of
    `spielman_wr`.
    """
    params = read_phydms_model_params(param)
    e_pw = [[] for x in range(3)]
    for i in range(3):
        params["phi{0}T".format(i)] = 1 - sum([params[x] for x in params.keys()
                                               if x.startswith("phi{0}"
                                               .format(i))])
        e_pw[i] = [params["phi{0}{1}".format(i, INDEX_TO_NT[x])] for x in
                   range(N_NT)]
    e_pw = scipy.array(e_pw)
    return phydmslib.models.YNGKP_M0(e_pw, 1, kappa=params["kappa"],
                                     omega=spielman_wr * params["omega"],
                                     mu=0.3, freeparams=['mu'])


# Calculations
def f_calculation(x, pr, Mt, model_name):
    """
    This function calcuates the value of `f` (the y-axis) for some amino-acid
    `x`, site `r`, and time `t`.

    The inputs for this function are the stationary state subsetted to site `r`
    and the transition matrix (M) subsetted to site `r` and time `t`.

    If the model is `YNGKP_M5`, then the input is a list of stationary states
    and transition matrices the length of the number of categories. The `f`
    value will be the average `f` value for the categories.

    >>> print(f_calculation(8, scipy.array([0.1,0.2,0.3,0.4]), \
    scipy.array([[1,-100,3,-100],[-100,-100,-100,-100],[9,-100,11,-100],\
    [-100,-100,-100,-100]]), "ExpCM"))
    6.4
    """
    target_codons = scipy.where(CODON_TO_AA == x)[0]  # amino-acid `x` codons
    if model_name == "GY94 + Gr":  # need to avg `f` over categories `k`
        prx = [pr[k][target_codons] for k in range(len(pr))]  # ss
        # transition to syn. codons
        Mtxy = [Mt[k][target_codons][:, target_codons] for k in range(len(Mt))]
        # sum the transition to syn codons
        Mtxy = [scipy.sum(Mtxy[k], axis=1) for k in range(len(Mt))]
        # multiply sum of transition to syn codons by ss and sum for each k
        f_rtx = [(Mtxy[k] * prx[k]).sum() for k in range(len(Mt))]
        return sum(f_rtx)/len(Mt)  # return the avg `f` for the categories `k`
    elif model_name in ["ExpCM", "GY94", "GY94 + wr"]:
        prx = pr[target_codons]  # stationary state
        Mtxy = Mt[target_codons][:, target_codons]  # transition to syn codons
        Mtxy = scipy.sum(Mtxy, axis=1)  # sum of transition to syn codons
        # multiply sum of the transition to syn codons by the ss and sum
        f_rtx = (prx * Mtxy).sum()
    else:
        raise ValueError("Cannot handle model {0}".format(model_name))
    return f_rtx


def get_pr(r, model, model_name):
    """
    The purpose of this function is to return the stationary state subsetted to
    site `r` correctly for each model.

    It returns a (61,61) array for every model except `YNGKP_M5` which it
    returns [(61,61) for k in n.cats].
    """
    if model_name == "GY94 + Gr":
        return [model.stationarystate(k)[r] for k in range(model.ncats)]
    elif model_name in ["ExpCM", "GY94", "GY94 + wr"]:
        return model.stationarystate[r]
    else:
        raise ValueError("Cannot handel model {0}".format(model_name))


def get_Mrt(r, t, model, model_name):
    """
    The purpose of this function is to return the transition matrix subsetted
    to site `r` and time `t`correctly for each model.

    This function also uses model.branchScale to correct `t` for each model.

    It returns a (61,61) array for every model except `YNGKP_M5` which it
    returns [(61,61) for k in n.cats].
    """
    if model_name == "GY94 + Gr":
        return [model.M(k, float(t/model.branchScale))[r] for k in
                range(model.ncats)]
    elif model_name in ["ExpCM", "GY94", "GY94 + wr"]:
        return model.M(float(t/model.branchScale))[r]
    else:
        raise ValueError("Cannot handel model {0}".format(model_name))


def main():
    """
    The basic workflow of this script is
    1. define the input files for the various models
    2. Create `ExCM` model
    3. Create `YNGKP_M0` and `YNGKP_M5` w/ the same # of sites as `ExpCM`
    4. Calculate the `spielman_wr` values from the `ExpCM`
    5. Loop through the sites and models. Create a `YNGKP_wr` for each site.
            Loop through times and amino-acids to calculate `f`.
    6. Output the dataframe.
    """

    # Set up the parameter files
    phydms_dir = "../HA/branch_lengths/phydms/"
    prefs_dir = "../HA/data/references/"
    ExpCM_modelparams_fname = ("{0}hybrid_lowH1_0_ExpCM_HA_hybridDoud_prefs_modelparams"
                               ".txt".format(phydms_dir))
    YNGKP_M0_modelparams_fname = ("{0}hybrid_lowH1_0_YNGKP_M0_modelparams"
                                  ".txt".format(phydms_dir))
    YNGKP_M5_modelparams_fname = ("{0}hybrid_lowH1_0_YNGKP_M5_modelparams"
                                  ".txt".format(phydms_dir))
    prefs_fname = "{0}HA_hybridDoud_prefs.csv".format(prefs_dir)
    if not os.path.isdir("outputs"):
        os.makedirs("outputs")
    # define the maximum amount of time
    max_time = 60

    # setup up the final dataframe
    df = {"Model": [], "Time": [], "f": [], "Site": []}

    # Make the models
    model_list = ["ExpCM", "GY94", "GY94 + wr", "GY94 + Gr"]
    models = {}
    if "ExpCM" in model_list:
        models["ExpCM"] = create_model_ExpCM(ExpCM_modelparams_fname,
                                             prefs_fname)
    else:
        raise ValueError("Must include `ExpCM`.")
    if "GY94" in model_list:
        models["GY94"] = create_model_YNGKP_M0(YNGKP_M0_modelparams_fname,
                                               models["ExpCM"].nsites)
    if "GY94 + Gr" in model_list:
        models["GY94 + Gr"] = create_model_YNGKP_M5(YNGKP_M5_modelparams_fname,
                                                    models["ExpCM"].nsites)
    if "GY94 + wr" in model_list:
        # calculate the spielman_wr values
        spielman_wr = models["ExpCM"].spielman_wr()
        wr = pd.DataFrame({"site": [x+1 for x in range(len(spielman_wr))],
                           "wr": spielman_wr})
        wr.to_csv("outputs/spielman_wr.csv", index=False)

    # Perform the calculations
    for r in range(models["ExpCM"].nsites):
        if r % 25 == 0:
            print(r)
        for model_name in model_list:
            if model_name != "GY94 + wr":
                model = models[model_name]
                _r = r  # dummy site to accomadate `YNGKP + wr`
            else:
                model = create_model_YNGKP_wr(YNGKP_M0_modelparams_fname,
                                              spielman_wr[r])
                _r = 0
            pr = get_pr(_r, model, model_name)
            for t in range(max_time):
                if t == 0:
                    f_rt = 1
                else:
                    f_rt = 0
                    Mrt = get_Mrt(_r, t, model, model_name)
                    for x in range(N_AA):
                        f_rt += f_calculation(x, pr, Mrt, model_name)
                df["Model"].append(model_name)
                df["Time"].append(t)
                df["f"].append(f_rt)
                df["Site"].append(r+1)
    df = pd.DataFrame(df)
    df.to_csv("outputs/expected_identity_given_time_t.csv", index=False)
    print("done.")


if __name__ == '__main__':
    main()
    import doctest
    doctest.testmod()
