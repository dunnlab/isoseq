# -*- coding: utf-8 -*-
import argparse
import ete3
import os

# from src.treeinform_collapse import load_trees
import matplotlib.pyplot as plt


def load_trees(gene_trees_folder):
    """
    Load tree paths from folder with gene trees.
    Returns a list with paths to each gene tree.

    Args:
    -- gene_trees_folder: path to folder with gene trees in newick format.

    """

    trees = []
    for file in os.listdir(gene_trees_folder):
        if file.endswith(".txt"):
            tree_path = os.path.join(gene_trees_folder, file)
            trees.append(tree_path)
    assert len(trees) != 0
    return trees


def branch_length_histogram(newicks):
    """
    # Agalma - an automated phylogenomics workflow
    # Copyright (c) 2012-2017 Brown University. All rights reserved.
    # Modified by Natasha Picciani on Oct 18

    Distribution of subtree branch lengths. Returns a dictionary
    containing subtree lengths and corresponding counts.

    Args:
    -- newicks: list containing the absolute paths to tree files
    """
    hist = {}
    for newick in newicks:
        tree = ete3.Tree(newick)
        outgroup = tree.get_midpoint_outgroup()
        if not outgroup is None:
            tree.set_outgroup(outgroup)
        # tree.convert_to_ultrametric(
        #     tree_length=1
        # ) # converts tree to ultrametric with length 1
        for node in tree.traverse(strategy="postorder"):
            if node.is_leaf():
                node.add_feature("branchlength", 0)
                node.add_feature("under", True)
            if not node.is_leaf():
                children = node.get_children()
                branchlength = (
                    children[0].get_distance(children[1])
                    + children[0].branchlength
                    + children[1].branchlength
                )
                node.add_feature("branchlength", branchlength)
                hist[branchlength] = hist.get(branchlength, 0) + 1
    return hist


def plot_histogram(histogram, binsize, threshold_value, outdir):
    """ "
    Plot distribution of subtree lengths.
    Returns png image with a histogram of the distribution.

    Args:
    -- histogram: dictionary containing subtree lengths and corresponding counts.
    -- outdir: output directory for image file named 'branch.length.hist.png'.
    """

    if threshold_value.is_integer():
        threshold_value_int = int(threshold_value)
        figname = f"{outdir}/branch.length.hist_{threshold_value_int}.png"
    else:
        figname = f"{outdir}/branch.length.hist_{threshold_value}.png"
    plt.hist(histogram, bins=binsize, edgecolor="k")
    plt.xlabel("Branch Length")
    plt.ylabel("Frequency")
    plt.title("Distribution of Subtree Branch Lengths")
    plt.axis([0, 30, 0, 20000])
    plt.axvline(threshold_value, color="red")
    plt.savefig(figname, facecolor="white")


def main(args):

    gene_trees_folder = args.gt
    outdir_name = args.o
    threshold_value = float(args.t)
    binsize = int(args.b)

    trees = load_trees(gene_trees_folder)
    histogram = branch_length_histogram(trees)
    plot_histogram(histogram, binsize, threshold_value, outdir_name)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Plot distribution of subtree branch lengths from several gene trees"
    )

    parser.add_argument(
        "-gt",
        "-genetrees",
        type=str,
        required=True,
        help="path to folder with gene trees produced with Orthofinder",
    )
    parser.add_argument(
        "-t",
        "-threshold",
        type=float,
        required=True,
        default="0.05",
        help="threshold of subtree length for collapsing gene variants",
    )
    parser.add_argument(
        "-b",
        "-binsize",
        type=int,
        required=True,
        default="10000",
        help="bin size for histogram",
    )

    parser.add_argument(
        "-o", "-outdir", type=str, required=False, default=".", help="output directory"
    )

    args = parser.parse_args()
    main(args)
