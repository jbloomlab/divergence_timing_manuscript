"""
The purpose of this script is to make the hybrid preference set between Mike
and Juhye's preferences.

SKH 20171219
"""

import pandas as pd
import numpy as np

def average_df(df1, df2):
    """
    >>> pd.DataFrame({"A":[0.45,1.1,3.5],"B":[-.5,1,4]}).equals(average_df(pd.DataFrame({"A":[1,2,3],"B":[-2,0,5]}), pd.DataFrame({"A":[-0.1,0.2,4],"B":[1,2,3]})))
    True
    """
    p = pd.Panel({n: df for n, df in enumerate([df1, df2])})
    average = p.mean(axis=0)
    return average

def main():
    # set up the file names
    doud_fname = "HA_Doud_prefs.csv"
    lee_fname = "HA_Lee_prefs.csv"
    hybrid_numbering = "hybrid_numbering.csv"
    outname = "HA_average_prefs.csv"

    # read in the files
    doud = pd.read_csv(doud_fname)
    lee = pd.read_csv(lee_fname)
    hybrid_numbers = pd.read_csv(hybrid_numbering)

    # subset to only the shared sites
    doud = doud[doud["site"].isin(hybrid_numbers["WSN"].tolist())]
    lee = lee[lee["site"].isin(hybrid_numbers["Perth"].tolist())]
    assert len(doud) == len(lee)

    # get rid of the sites column
    doud = doud.drop("site", axis=1)
    lee = lee.drop("site", axis=1)

    # average the values
    average = average_df(doud, lee)
    assert np.allclose(average.sum(axis=1), [1] * len(average))

    # add back in the site column
    average["site"] = [x+1 for x in range(len(average))]
    cols = ["site"] + [x for x in average.columns.values if x not in ["site", "index"]]
    average = average[cols]

    # output the file
    average.to_csv(outname, index=False)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    main()
