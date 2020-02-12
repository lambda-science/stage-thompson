# Script permettant de co-localiser les petit exons/introns avec la zone génomique/CDS du mismatch
# INPUT
# arg1: mismatch_exon_file: fichier contenant les sequences des exons localisés au niveau des mismatch
# arg2: mismatch_intron_file: pareil mais pour les introns
# arg3: chemin+nom du fichier d'output pour la seq de CDS mismatch
# arg4: chemin+nom du fichier d'output pour la seq de genomique mismatch
# OUTPUT
# Deux fichiers arg3 et arg4

import pandas as pd
import sys


def fasta2List(pathFasta):
    # Importer fichier fasta en dictionnaire
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


def getCDSMismatchSeq(mismatch_exon_file, CDS_seq_out_file):
    # Récup CDS et genomic que au niveau du mismatch
    # Bug to fix: problême duplication dans fichier d'erreur = CDS deux fois plus longue que prévu :(
    myExonicMap = pd.read_csv(mismatch_exon_file, sep="\t", header=None)
    UniprotID = set(myExonicMap.iloc[:, 0].to_list())
    UniprotID = list(UniprotID)

    with open(CDS_seq_out_file, "w") as CDS_file:
        for i in UniprotID:
            subset = myExonicMap.loc[myExonicMap[0] == i]
            myCDS = ''.join(subset[6].to_list())
            CDS_file.write(">"+subset.iloc[0, 1] +
                           " "+i+" CDS_mismatch"+"\n"+myCDS+"\n")


def getGenomicMismatchSeq(mismatch_exon_file, mismatch_intron_file, genomic_seq_out_file):
    # Récup genomic que au niveau du mismatch
    exon_file = pd.read_csv(mismatch_exon_file, sep="\t", header=None)
    intron_file = pd.read_csv(mismatch_intron_file, sep="\t", header=None)
    UniprotID = set(exon_file.iloc[:, 0].to_list())
    UniprotID = list(UniprotID)

    with open(genomic_seq_out_file, "w") as CDS_file:
        for i in UniprotID[:]:
            my_CDS = []
            subset_exon = exon_file.loc[exon_file[0] == i]
            subset_intron = intron_file.loc[intron_file[0] == i]
            for n in range(0, len(subset_exon.index)):
                my_CDS.append(subset_exon.iloc[n, 6])
                try:
                    my_CDS.append(subset_intron.iloc[n, 6])
                except:
                    pass
            CDS_file.write(
                ">"+subset_exon.iloc[0, 1]+" "+i+" Genomic_mismatch"+"\n"+''.join(my_CDS)+"\n")


if __name__ == "__main__":
    mismatch_exon_file = sys.argv[1]
    mismatch_intron_file = sys.argv[2]
    CDS_out_file = sys.argv[3]
    genomic_out_file = sys.argv[4]

    getCDSMismatchSeq(mismatch_exon_file, CDS_out_file)
    getGenomicMismatchSeq(mismatch_exon_file,
                          mismatch_intron_file, genomic_out_file)
