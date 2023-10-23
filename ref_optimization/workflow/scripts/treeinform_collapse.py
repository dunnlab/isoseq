# -*- coding: utf-8 -*-
import argparse
import ete3
import os
import operator
from collections import defaultdict
from Bio import SeqIO
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


def identify_candidate_variants(newicks, threshold, species_of_interest):
    """
    # Agalma - an automated phylogenomics workflow
    # Copyright (c) 2012-2017 Brown University. All rights reserved.
    # Modified by Natasha Picciani on Oct 18

    For each gene tree, identify candidates for variants
    of the same gene based on a threshold for branch lengths of subtrees.
    Return a list containing sets of variants of the same gene
    for the species of interest.

    Obs: transcript id must follow pattern "species_pep_transcript".
    Otherwise, change separator in line 66

    Args:
    -- newicks: list containing the absolute paths to tree files.
    -- threshold: branch length threshold under which sequences will be
    considered variants of the same gene
    -- species_of_interest: species considered for retrieving gene
    variants. Species name format must match that in gene trees;
    e.g. Podocoryna_carnea.
    """

    candidates_specific = []
    for newick in newicks:
        tree = ete3.Tree(newick)
        outgroup = tree.get_midpoint_outgroup()
        if not outgroup is None:
            tree.set_outgroup(outgroup)
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
                if branchlength < threshold:  # adds flag for if node.bl under threshold
                    node.add_feature("under", True)
                else:
                    node.add_feature("under", False)

        # yield candidates for subtrees with branch length < threshold
        for node in tree.traverse(
            strategy="levelorder", is_leaf_fn=operator.attrgetter("under")
        ):
            if node.branchlength != 0 and node.under == True:
                candidates = defaultdict(
                    set
                )  # default factory is a set so no duplicate values are added
                for leaf in node.get_leaves():
                    species, _, model_id = leaf.name.partition(
                        "_pep_"
                    )  # model id is gene model id
                    model_id = str(model_id)
                    candidates[species].add(model_id)
                for child in node.get_children():
                    child.checked = True
                if (
                    len(candidates[species_of_interest]) > 1
                ):  # get variants from species of interest
                    # yield (candidates[species_of_interest])
                    candidates_specific.append(candidates[species_of_interest])
    return candidates_specific


def create_protein_length_dictionary(peptides):
    """
    Returns a dictionary containing protein sequence records as keys
    and their corresponding lengths as value.

    Args:
    -- peptides: path to fasta file containing protein sequences
    """
    protein_sequences = {}
    protein_lengths = {}

    with open(peptides) as input_seqs:
        protein_sequences = SeqIO.to_dict(SeqIO.parse(input_seqs, "fasta"))
        for protein in protein_sequences:
            protein_lengths[protein] = len(protein_sequences[protein])
    return protein_lengths


def fix_special_characters(protein_lengths):
    """
    Replaces special characters ':', '(', and ')' with '_' from dictionary with
    sequence ids loaded from protein fasta file generated with transdecoder.
    Returns a dictionary with fixed protein names matching those from
    orthofinder gene tree tips.

    Args:
    -- protein_lengths: dictonary with protein sequences and lengths
    generated with create_protein_length_dictonary.

    """
    new_protein_lengths = (
        {}
    )  # make names match those from orthofinder tips where some special characters were replaced with _
    for key, value in protein_lengths.items():
        new_key = key.replace(":", "_").replace("(", "_").replace(")", "_")
        new_protein_lengths[new_key] = value
    return new_protein_lengths


def retrieve_longest_variant(candidate_variants, new_protein_lengths):
    """
    Retrieve the longest variant of a gene. If several sequences have equal lengths,
    the longest is randomly chosen. Returns a list with names of selected variants
    across all sets of variants.

    Args:
    -- candidate_variants: list from identify_candidate_variants
    -- new_protein_lengths: dictionary of protein lengths generated after running fix_special_characters
    """

    selected_variants = []
    for set_of_variants in candidate_variants:
        variantset = {}
        for sequence in set_of_variants:
            length = new_protein_lengths[sequence]
            variantset[sequence] = length
        selected_variants.append(max(variantset, key=variantset.get))
    return selected_variants


def filter_sequences(candidate_variants, selected_variants, peptides):
    """
    Generate a list of sequences one transcript variant, the longest, per gene.

    Args:
    -- candidate_variants: list from identify_candidate_variants
    -- selected_variants: list of longest sequence per set of variants from retrieve_longest_variant
    -- peptides: path to fasta file with protein sequences
    """
    all_variants = []
    for set_of_variants in candidate_variants:
        for sequence in set_of_variants:
            all_variants.append(sequence)

    exclude = []
    for variant in all_variants:
        if variant not in selected_variants:
            exclude.append(variant)

    filtered_sequences = []
    with open(peptides) as input_seqs:
        for record in SeqIO.parse(input_seqs, "fasta"):
            record.id = record.id.replace(":", "_").replace("(", "_").replace(")", "_")
            if record.id in exclude:
                continue
            else:
                filtered_sequences.append(record)
    return filtered_sequences


def export_variants(candidate_variants, species_of_interest, outdir):
    """
    Export csv file with one set of variants per line
    """
    outfile = f"{outdir}/{species_of_interest}.variants.csv"
    with open(outfile, "w") as out:
        for set_of_variants in candidate_variants:
            lst = ",".join(
                str(s) for s in set_of_variants
            )  # get rid of set curly brackets and use "," as separator
            out.write(lst + "\n")


def main(args):

    protein_file = args.s
    gene_trees_folder = args.gt
    outdir_name = args.o
    threshold_value = float(args.t)
    species_of_interest = args.sp

    trees = load_trees(gene_trees_folder)
    candidate_variants = identify_candidate_variants(
        trees, threshold_value, species_of_interest
    )
    export_variants(candidate_variants, species_of_interest, outdir_name)
    new_protein_lengths = fix_special_characters(
        create_protein_length_dictionary(protein_file)
    )
    selected_variants = retrieve_longest_variant(
        candidate_variants, new_protein_lengths
    )
    final_sequences = filter_sequences(
        candidate_variants, selected_variants, protein_file
    )

    # export fasta file after filtering sequences
    outfile2 = f"{outdir_name}/{species_of_interest}.collapsed.fasta"
    with open(outfile2, "w") as out2:
        for sequence in final_sequences:
            SeqIO.write(sequence, out2, "fasta")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Collapse protein fasta file by selecting one sequence variant per gene based on the topology of gene trees"
    )

    parser.add_argument(
        "-s",
        "-sequence",
        type=str,
        required=True,
        help="protein fasta file from species of interest used when generating gene trees",
    )
    parser.add_argument(
        "-gt",
        "-genetrees",
        type=str,
        required=True,
        help="path to folder with gene trees produced with Orthofinder",
    )
    parser.add_argument(
        "-sp",
        "-species",
        type=str,
        required=True,
        default=".",
        help="species for which to retrieve gene variants",
    )
    parser.add_argument(
        "-t",
        "-threshold",
        type=float,
        required=True,
        default="0.05",
        help="subtree length value under which sequences are considered variants",
    )
    parser.add_argument(
        "-o", "-outdir", type=str, required=False, default=".", help="output directory"
    )

    args = parser.parse_args()
    main(args)
