# Script permettant de tenter de corriger un mismatch en alignant le peptide human a la seq génomic et en cherchant le meilleur match
# Input:
#     arg1 Error_file fichier contenant toute les erreurs de type mismatch (script julie)
#     arg2 my_Genomic fichier fasta contenant les séquences génomiques des transcripts
#     arg3 out_folder: dossier de sortie qui va contenir les fichier fasta et mafft d'ailgnement
#     arg4 mafft_path: chemin vers le programme mafft
#     arg5 results_file: chemin vers le fichier final contenant les matchs
# Output
# Dossier fasta/ mafft/ et fichier results_file

import json
import pandas as pd
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
import os


def fasta2List(pathFasta):
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


def translateAndAlign(Error_file, out_folder, mafft_path):

    # Pour chaque mismatch du fichier d'erreur: on cherche la séquence humaine en face du mismatch. On prend le séquence génomique de la protéine du primate
    # concerné et on la traduit dans les trois cadre de lecture. On cherche si la séquence humaine peut être trouvée dans une des séquence génomique traduite
    # Si oui on écrit dans un fichier  les informations du mismatch et la position génomique permetttant de le régler.
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
        os.mkdir(out_folder+"/fasta/")
        os.mkdir(out_folder+"/mafft/")

    for index, row in Error_file.iloc[:, :].iterrows():
        fasta_name = row[0][20:-6]
        prot_name = row[2]
        try:
            # For mismatch: col 5 and 6 are human peptide, 3 and 4 are primate
            human_start = row[5]
            human_stop = row[6]
        except:
            # For Type 2 error col 3 and 4 are human peptide, primate are absent
            human_start = row[3]
            human_stop = row[4]

        Prot_list = fasta2List("../data/raw/uniprot-sequence/"+fasta_name)
        prot_HumanRef = [val for key, val in Prot_list.items()
                         if row[0][20:-15] in key]
        genomic_Seq = [val for key, val in my_Genomic.items()
                       if prot_name in key]
        if genomic_Seq == []:
            pass

        elif prot_HumanRef == []:
            pass
        else:
            peptide_Ref = (row[0][20:-15], prot_HumanRef[0]
                           [human_start:human_stop])

            # Translation & alignement part
            for i in range(3):
                f = open(out_folder+"fasta/"+str(index) +
                         "_"+peptide_Ref[0]+"_"+str(i)+".fasta", "w")
                genomic_prot = Seq(genomic_Seq[0][i:], generic_dna)
                f.write(">"+prot_name+" genomic frame "+str(i) +
                        "\n"+str(genomic_prot.translate())+"\n")
                f.write(
                    ">"+peptide_Ref[0]+" Human Peptide query"+"\n"+str(peptide_Ref[1])+"\n")
                f.close()

                # Alignement MAFFT avec paramètre strictes sur l'ouverte et l'extension de gap
                mafft_cline = MafftCommandline(mafft_path,
                                               input=out_folder+"/fasta/" +
                                               str(index)+"_" +
                                               peptide_Ref[0]+"_" +
                                               str(i)+".fasta",
                                               op=5.0, ep=10.0)
                mafft_results = mafft_cline()
                with open(out_folder+"mafft/"+str(index)+"_"+peptide_Ref[0]+"_"+str(i)+".fasta.mafft", "w") as handle:
                    handle.write(mafft_results[0])


def selectMatchInAlignement(mafft_align_folder, results_file):
    f = open(results_file, "w")
    align_files = os.listdir(mafft_align_folder)
    # Scan de la séquence aligné peptide humain pour trouver la position de début et de fin sur l'alignement
    f.write(
        "Match\tPrimate\tHuman\tStartPos\tStopPos\tSequence_Primate\tSequence_Humaine\n")
    for align_file in align_files:
        align = AlignIO.read(mafft_align_folder+align_file, "fasta")
        start = 0
        last_amino_index = 0  # Debug ?
        found_Start = False
        for index, amino in enumerate(align[1, :].seq):
            if amino.isalpha() and found_Start == False:
                start = index
                found_Start = True
            if found_Start and amino.isalpha():
                last_amino_index = index+1

        peptide_genomique = align[0, start:last_amino_index].seq
        peptide_humain_query = align[1, start:last_amino_index].seq

        # Boucle qui scan l'identité entre la paire de peptide
        identity_count = 0
        gap_count = 0
        for i in range(len(peptide_humain_query)):
            if peptide_genomique[i] == peptide_humain_query[i]:
                identity_count += 1
            if (peptide_genomique[i] == "-") or (peptide_humain_query[i] == "-"):
                gap_count += 1
        identity_percentage = identity_count/len(peptide_humain_query)

        if gap_count < 10 and identity_percentage > 0.8:
            f.write("Match\t"+align[0].id+"\t"+align[1].id+"\t"+str(start*3)+"\t"+str(last_amino_index*3)+"\t"
                    + str(peptide_genomique)+"\t"+str(peptide_humain_query)+"\n")
        else:
            pass
    f.close()


if __name__ == "__main__":
    # Importation de toutes les données utilisées
    Error_file = pd.read_csv(sys.argv[1], sep=" ", header=None)
    my_Genomic = fasta2List(sys.argv[2])
    out_folder = sys.argv[3]
    mafft_align_folder = out_folder+"mafft/"
    mafft_path = sys.argv[4]
    results_file = sys.argv[5]
    if len(os.listdir(out_folder+"mafft")) < 60000:  # DEBUGGIN CONDITION
        translateAndAlign(Error_file, out_folder, mafft_path)
    else:
        selectMatchInAlignement(mafft_align_folder, results_file)
