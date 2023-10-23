import pandas as pd

targets = (
    pd.read_csv(config["targets"], sep="\t", dtype={"species": str})
    .set_index("species", drop=False)
    .sort_index()
)

ensembl_targets = targets.loc[lambda targets: targets['source'] == "ensembl"]
ensemblgenomes_targets = targets.loc[lambda targets: targets['source'] == "ensemblgenomes"]
other_targets = targets.loc[lambda targets: targets['source'] == "other"]
other_gz_targets = targets.loc[lambda targets: targets['source'] == "other_gz"]
gdrive_targets = targets.loc[lambda targets: targets['source'] == "gdrive"]

def get_sequence(wildcards, type):
    """
    Return path for downloading sequence file.
    Wildcards: target species
    Type: defined in file "download_targets.tsv" (ensembl, ensembl genomes, other, other gz or gdrive)
    """
    species_units = targets.loc[wildcards.species]
    if type == "ensembl":
        if    species_units["source"] == "ensembl":
            return species_units["file"]
    if type == "ensemblgenomes":
        if    species_units["source"] == "ensemblgenomes":
            return species_units["file"]
    if type == "other":
        if species_units["source"] == "other":
            return species_units["file"]
    if type == "other_gz":
        if species_units["source"] == "other_gz":
            return species_units["file"]
    if type == "gdrive":
        if species_units["source"] == "gdrive":
            return species_units["file"]

rule get_ensembl:
    """
    Download species protein file from ENSEMBL.
    """
    output:
        "resources/sequences/ensembl/{species}.pep.fasta.gz" # this renames the original file and keep only the species names as given in the species column of targets dataframe
    params:
        lambda wc: get_sequence(wc, type="ensembl")
    shell:
        "wget http://ftp.ensembl.org/{params} -O {output}"

rule get_ensemblgenomes:
    """
    Download species protein file from ENSEMBL Genomes.
    """
    output:
        "resources/sequences/ensemblgenomes/{species}.pep.fasta.gz" # this renames the original file and keep only the species names as given in the species column of targets dataframe
    params:
        lambda wc: get_sequence(wc, type="ensemblgenomes")
    shell:
        "wget http://ftp.ensemblgenomes.org/{params} -O {output}"

rule get_other:
    """
    Download species protein file from random websites.
    """
    output:
        "resources/sequences/other/{species}.pep.fasta"
    params:
        lambda wc: get_sequence(wc, type="other")
    shell:
        """
        wget --no-check-certificate {params} -O {output}
        """

rule get_other_gz:
    """
    Download compressed (gzipped) species protein file from random websites.
    """
    output:
        "resources/sequences/other_gz/{species}.pep.fasta.gz"
    params:
        lambda wc: get_sequence(wc, type="other_gz")
    shell:
        """
        wget --no-check-certificate {params} -O {output}
        """

rule get_gdrive:
    """
    Download species protein file from Google Drive.
    """
    output:
        "resources/sequences/gdrive/{species}.pep.fasta"
    params:
        lambda wc: get_sequence(wc, type="gdrive")
    shell:
        """
        wget -O {output} {params}
        """
#filename and path is embedded in the gdrive command listed in the file column of target tsv
