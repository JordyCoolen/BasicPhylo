import dendropy
import argparse

######
INFO = "Create dendogram from distance matrix"
__author__ = "J.P.M. Coolen"
__version__ = 0.1
######

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="full path of input distance matrix"),
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Prefix of output")
    parser.add_argument("--outgroup", type=str, required=False,
                        help="select outgroup for tree", default="Hoyosella_altamirensis.fasta")
    parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s {version}".format(version=__version__))

    # parse all arguments
    args = parser.parse_args()

    return args

if __name__ == "__main__":
    # load arguments set global workdir
    args = parse_args()

    # Load the distance matrix from a file
    pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
            args.input,
            is_first_row_column_names=True,
            is_first_column_row_names=True,
            is_allow_new_taxa=True,
            delimiter="\t",
            )

    # Calculate UPGMA from distance matrix
    upgma_tree = pdm.upgma_tree()

    # select outgroup and reroot tree
    try:
        outgroup_node = upgma_tree.find_node_with_taxon_label(args.outgroup)
        upgma_tree.to_outgroup_position(outgroup_node, update_bipartitions=False)
    except AttributeError:
        print('Wrong outbreak, please correct. Rerooting will be skipped')

    # Print data and tree (for debugging)
    #print(upgma_tree.as_string("nexus"))
    #print(upgma_tree.as_ascii_plot())

    # write newick tree to file
    upgma_tree.write_to_path(f"{args.output}.tre", schema="newick")

    print("Finished")