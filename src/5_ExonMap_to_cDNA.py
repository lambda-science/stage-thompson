# # Script permettant de générer les cDNA à partir de la carte exonique des gènes
# Input: Exon_map.tab
# Output: cdna_seq.fa

import pandas as pd
import sys


def writecDNASeq(myExonicMap, cDNAFile):
    # Function: write a fasta file containing cDNA sequence for each transcript of interest.
    # Parameters:
    # 		myExonicMap: (str) path to the exonic map file to read
    #       cDNAFile: (str) path to the cDNA fasta file to write
    # Return: None. Write a file at cDNAFile path.
    UniprotID = set(myExonicMap.iloc[:, 0].to_list())
    UniprotID = list(UniprotID)
    with open(mycDNAFile, "w") as cdna_file:
        for i in UniprotID:
            subset = myExonicMap.loc[myExonicMap[0] == i]
            myCDS = ''.join(subset[6].to_list())
            cdna_file.write(
                ">"+subset.iloc[0, 1]+" "+i+" cDNA"+"\n"+myCDS+"\n")


if __name__ == "__main__":
    myExonicMap = pd.read_csv(sys.argv[1], sep="\t", header=None)
    mycDNAFile = sys.argv[2]
    writecDNASeq(myExonicMap, mycDNAFile)
