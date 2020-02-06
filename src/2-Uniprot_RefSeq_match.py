# # Script permettant de chercher l'existance de séquences à mismatch Uniprot dans la base locale RefSeq
# Input:
#     fasta contenant les séquence Uniprot d'intérêt  (my_Query)
#     fasta contenant toute les protéines présente dans RefSeq (my_DB)
# Output:  fichier match.txt une ligne par match entre les deux fichiers

import sys


def fasta2List(pathFasta):
    # Lis un fasta et retourne un dictionnaire de type: sequence:titre
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
    dictionary = dict(zip(seq, title))
    return dictionary


def importQuery(pathFasta):
    # Lis un fasta et retourne un dictionnaire de type: titre:sequence
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


if __name__ == "__main__":
    my_DB = fasta2List(sys.argv[1])
    my_Query = importQuery(sys.argv[2])

    f = open("raw/uniprot-error-mismatch/uniprot_refseq_match.out", "w")
    for i, j in my_Query:
        try:
            match = my_DB[j]
            f.write("Match found: " + str(my_DB[j])+" "+str(i))
        except:
            pass
