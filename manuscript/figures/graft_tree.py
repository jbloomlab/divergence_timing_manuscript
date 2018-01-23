"""
The purpose of this script is to graft a small tree ("cion") onto a larger
tree ("stock").

SKH 20180120
"""

from ete3 import Tree
import argparse
import sys


class ArgumentParserNoArgHelp(argparse.ArgumentParser):
    """Like *argparse.ArgumentParser*, but prints help when no arguments."""
    def error(self, message):
        sys.stderr.write('error: {0}\n\n'.format(message))
        self.print_help()
        sys.exit(2)


def ParseArguments():
    """Argument parser for script."""
    parser = ArgumentParserNoArgHelp(
            description='Graft trees',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            )
    parser.add_argument('stock', help='File name of stock tree.')
    parser.add_argument('cion', help='File name of cion tree.')
    parser.add_argument('out_prefix', help='Prefix for the final tree file')
    return parser


def tree_graft(stock, cion):
    names = []
    for node in cion.traverse("postorder"):
        if len(node.name) >= 1:
            names.append(node.name)
    target_node = stock.get_common_ancestor(names)
    for child in target_node.get_children():
        child.detach()
    cion_root = cion.get_tree_root()
    for child in cion_root.get_children():
        target_node.add_child(child)
    return stock


def main():
    args = vars(ParseArguments().parse_args())
    stock = args["stock"]
    cion = args["cion"]

    catch = tree_graft(Tree(stock), Tree(cion))
    catch.write(format=1, outfile="{0}.newick".format(args["out_prefix"]))


if __name__ == '__main__':
    main()
