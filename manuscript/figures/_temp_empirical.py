"""
The purpose of this script is to manipulate tree figures to make the
alternative tree figures.

SKH 20180131
"""

from ete3 import Tree
import glob
import os

def main():
    tree_list = glob.glob("../../analysis/HA/branch_lengths/phydms/hybrid_high_0_*_tree.newick")
    tree_list = [x for x in tree_list if "randomized" not in x]
    tree_list = [x for x in tree_list if "averaged" not in x]
    tree_list = [x for x in tree_list if "RAxML" not in x]
    assert len(tree_list) == 8
    for tree in tree_list:
        outName = "_temp_{0}_{1}.newick".format(os.path.basename(tree).split(".")[0], "out")
        t = Tree(tree)
        R = t.get_midpoint_outgroup()
        t.set_outgroup(R)
        t.write(format=1, outfile=outName)

if __name__ == '__main__':
    main()
