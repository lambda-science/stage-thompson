# # Script permettant de tenter de corriger un mismatch en cherchant le peptide human dans la séquence génomique de la protéine mal prédite
# Input:
#     arg1 Error_file fichier contenant toute les erreurs de type mismatch (script julie)
#     arg2 my_Genomic fichier fasta contenant les séquences génomiques des transcripts
#     arg3 results_file fichier ou écrire l'output

# Importer les lib et fonction
import json
import pandas as pd
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


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


def translationCorrectMismtach(results_file, Error_File, my_Genomic):
    # Function: Try to correct a mismatch or deletion by translating the genomic sequence of the transcript and looking for the human peptid.
    # Parameters:
    # 		results_file: (str) path to the results file to write
    #       Error_file: (str) path to the error file to read
    #       my_Genomic: (dict) dictionnary containing genomique sequence of transcript of interest
    # Return: None. Write a file at results_file path.
    # Description:
    # For each mismatch of the error file we look for the human peptide creating the mismatch with the primate. We take the genomique sequence
    # of the mismatching primate protein and translate it in the three frames. We then search if the human sequence can be found in the primate
    # translated genomic one. If there is a match we write the informations in the results folder with the genomic position involved.
    Error_File = pd.read_csv(Error_File, sep=" ", header=None)
    ff = open(results_file, "w")
    ff.write("Match\tPrimate\tHuman\tStartPos\tStopPos\tSequence\n")
    for index, row in Error_file.iloc[:, :].iterrows():
        fasta_name = row[0][20:-6]
        prot_name = row[2]
        human_start = row[5]
        human_stop = row[6]

        Prot_list = fasta2List("uniprot/"+fasta_name)
        prot_HumanRef = [val for key, val in Prot_list.items()
                         if row[0][20:-15] in key]
        genomic_Seq = [val for key, val in my_Genomic.items()
                       if prot_name in key]

        if prot_HumanRef == []:
            pass
        else:
            peptide_Ref = (row[0][20:-15], prot_HumanRef[0]
                           [human_start:human_stop])
            # Translation part
            querySeq = Seq(peptide_Ref[1])
            for i in range(3):
                genomic_prot = Seq(genomic_Seq[0][i:], generic_dna)
                prot_trial = genomic_prot.translate()
                if (querySeq in prot_trial):
                    start = prot_trial.find(querySeq)*3
                    stop = start+len(querySeq)*3
                    ff.write("Match\t"+prot_name+"\t"+row[0][20:-15]+"\t"+str(
                        start)+"\t"+str(stop)+"\t"+str(querySeq)+"\n")
                    match_counter += 1

    ff.close()


if __name__ == '__main__':
    Error_file = sys.argv[1]
    my_Genomic = fasta2List(sys.argv[2])
    translationCorrectMismtach(sys.argv[3], Error_file, my_Genomic)
