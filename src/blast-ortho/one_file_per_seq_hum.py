# %%

import pandas as pd
from os import listdir
from os.path import isfile, join
import json

# %%


def fasta2List(pathFasta):
    # Function: Convert fastafile to dictionnary with title:seq structure
    # Parameters:
    # 		pathFasta: (str) path to the fasta file
    # Return:
    # 		dictionary: (dict) dictionnary title:seq structure
    f = open(pathFasta, "r")
    title = []
    seq = []
    seq_temp = []
    for line in f:
        if line[0] == ">":
            seq.append(''.join(seq_temp).replace("\n", ""))
            title.append(line.replace("\n", ""))
            seq_temp = []
        else:
            seq_temp.append(line)
    seq.append(''.join(seq_temp).replace("\n", ""))
    seq.pop(0)
    dictionary = dict(zip(title, seq))
    return dictionary


def write_hum_prot(fasta_dict):
    for i, j in fasta_dict.items():
        prot_name = i.split("|")[1]
        f = open("../../data/raw/uniprot-blast/" +
                 prot_name+".id", "w", newline='\n')
        f.write(str(i)+"\n")
        f.write(str(j)+"\n")
        f.close()

        g = open("../../data/raw/refseq-blast/" +
                 prot_name+".id", "w", newline='\n')
        g.write(str(i)+"\n")
        g.write(str(j)+"\n")
        g.close()

# %%


if __name__ == "__main__":
    fasta_dict = fasta2List("../../temp/UP000005640_9606.fasta")
    write_hum_prot(fasta_dict)
