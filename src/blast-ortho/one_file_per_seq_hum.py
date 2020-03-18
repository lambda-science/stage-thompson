import pandas as pd
from os import listdir
from os.path import isfile, join
import json


def write_orthologue_ID(gene_names):
    # Function: Separate row of orthologous uniprot ID of all primate table to one row per file containing each ID.
    # Parameters:
    # 		gene_names: (dataframe) dataframe containing uniprotID for each human prot per row and each primates by columns
    # Return: None. Write a file for each row in /raw/uniprot-sequence/
    # Description:
    for i in range(20265):
        f = open("../../data/raw/uniprot-blast/" +
                 str(gene_names.iloc[i, 0])+".id", "w", newline='\n')
        f.write(str(gene_names.iloc[i, 0])+"\n")
        f.close()

        g = open("../../data/raw/refseq-blast/" +
                 str(gene_names.iloc[i, 0])+".id", "w", newline='\n')
        g.write(str(gene_names.iloc[i, 0])+"\n")
        g.close()


if __name__ == "__main__":
    gene_names = pd.read_csv(
        "../../data/raw/orthoinspector-json-processed/allProteineName.txt", sep="\t")
    write_orthologue_ID(gene_names)
