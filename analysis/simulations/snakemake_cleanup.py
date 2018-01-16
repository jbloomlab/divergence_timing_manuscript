#!/python

"""
This script cleans up the out files from the Halern-Bruno-like simulations `Snakefile`.

This script is greedy and not robust.

SKH 20170420
"""

import glob
import os
import sys
import argparse

class ArgumentParserNoArgHelp(argparse.ArgumentParser):
    """Like *argparse.ArgumentParser*, but prints help when no arguments."""
    def error(self, message):
        sys.stderr.write('error: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


def ParseArguments():
    """Argument parser for script."""
    parser = ArgumentParserNoArgHelp(
            description='Clean up files from simulation',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            )
    parser.set_defaults(tree=False)
    parser.add_argument('--tree', dest='tree', action=\
            'store_true', help="Remove temp tree files.")
    parser.set_defaults(alignment=False)
    parser.add_argument('--alignment', dest='alignment', action=\
            'store_true', help="Remove alignment files.")
    parser.set_defaults(phydms=False)
    parser.add_argument('--phydms', dest='phydms', action=\
            'store_true', help="Remove phydms files.")
    parser.set_defaults(slurm=False)
    parser.add_argument('--slurm', dest='slurm', action=\
            'store_true', help="Remove slurm files.")
    return parser

def delete_files(keyword):
    for f in glob.glob("{0}*".format(keyword)):
        if os.path.isfile(f):
            os.remove(f)

def main():
    keyword_dictionary = {"tree":"_tempTree_", "alignment":"*simulatedalignment.fasta", "phydms":"phydms_", "slurm":"slurm-"}
    args = vars(ParseArguments().parse_args())
    print("Read the following command line arguments:")
    print("\t{0}".format("\n\t".join(["{0} = {1}".format(key, value) for (key, value) in args.items()])))

    if True in args.values():
        print("Cleaning up files from simulations")
        for arg in args.keys():
            if args[arg] == True:
                print("Removing {0}".format(arg))
                delete_files(keyword_dictionary[arg])
    print("done.")

if __name__ == '__main__':
 	main() # run the script
