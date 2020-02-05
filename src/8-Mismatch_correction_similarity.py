# Script permettant de tenter de corriger un mismatch en alignant le peptide human a la seq génomic et en cherchant le meilleur match
# Input:   
#     ID_file fichier de correspondance UniprotID - Ensembl Trasncript ID  
#     Error_file fichier contenant toute les erreurs de type mismatch (script julie)  
#     Exon_file fichier contenant les informations sur les positions des exons des transcripts  
#     my_Genomic fichier fasta contenant les séquences génomiques des transcripts  
#     my_cDNA fichier fasta contenant les séquences cDNA des transcripts  
#     my_CDS fichier fasta contenant les séquences CDS des transcripts 

# Importer les lib et fonction
import json
import pandas as pd
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

# Importation de toutes les données utilisées
ID_file = pd.read_csv("../raw/uniprot-exon-map/transcript_ensembl.tab", sep = "\t")
Error_file = pd.read_csv("../raw/uniprot-exon-map/uniprot_new_errors_filt.txt", sep=" ", header=None)
Exon_file = pd.read_csv("../raw/uniprot-exon-map/Exon_map.tab", sep="\t", header=None)
my_Genomic = fasta2List("../raw/uniprot-exon-map/genomics_new.fa")
my_cDNA = fasta2List("../raw/uniprot-exon-map/cdna_seq.fa") 
my_CDS = fasta2List("../raw/uniprot-exon-map/cds_new.fa")

# Pour chaque mismatch du fichier d'erreur: on cherche la séquence humaine en face du mismatch. On prend le séquence génomique de la protéine du primate
# concerné et on la traduit dans les trois cadre de lecture. On cherche si la séquence humaine peut être trouvée dans une des séquence génomique traduite
# Si oui on écrit dans un fichier  les informations du mismatch et la position génomique permetttant de le régler. 
for index, row in Error_file.iloc[:1,:].iterrows():
    fasta_name = row[0][20:-6]
    prot_name = row[2]
    error_start = row[3]
    error_stop = row[4]
    human_start = row[5]
    human_stop = row[6]

    Prot_list = fasta2List("../raw/uniprot-sequence/"+fasta_name)
    prot_HumanRef = [val for key, val in Prot_list.items() if row[0][20:-15] in key]
    genomic_Seq = [val for key, val in my_Genomic.items() if prot_name in key]
    
    if prot_HumanRef == []:
        pass
    else:
        peptide_Ref = (row[0][20:-15], prot_HumanRef[0][human_start:human_stop])
        
        # Translation & alignement part
        querySeq = Seq(peptide_Ref[1])
        prot_trial = []
        for i in range(3):
            f = open("../temp/"+str(index)+"_"+peptide_Ref[0]+"_"+str(i)+".fasta", "w")
            genomic_prot = Seq(genomic_Seq[0][i:], generic_dna)
            f.write(">"+prot_name+" genomic frame "+str(i)+"\n"+str(genomic_prot.translate())+"\n")
            f.write(">"+peptide_Ref[0]+" Human Peptide query"+"\n"+str(peptide_Ref[1])+"\n")
            f.close()
            
            # Alignement MAFFT avec paramètre strictes sur l'ouverte et l'extension de gap
            mafft_cline = MafftCommandline(input="../temp/"+str(index)+"_"+peptide_Ref[0]+"_"+str(i)+".fasta", op=5.0, ep=10.0)
            mafft_results = mafft_cline()
            with open("../temp/"+str(index)+"_"+peptide_Ref[0]+"_"+str(i)+".fasta.mafft", "w") as handle:
                handle.write(mafft_results[0])

from Bio import AlignIO
import os

f = open("../raw/correction-pairwise/translation_match.tab", "w")
align_files = os.listdir("../raw/correction-pairwise/mafft/")
# Scan de la séquence aligné peptide humain pour trouver la position de début et de fin sur l'alignement
f.write("Match\tPrimate\tHuman\tStartPos\tStopPos\tSequence_Primate\tSequence_Humaine\n")
for align_file in align_files:
    align = AlignIO.read("../raw/correction-pairwise/mafft/"+align_file, "fasta")
    start = 0
    found_Start = False
    for index, amino in enumerate(align[1, :].seq):
        if amino.isalpha() and found_Start == False:
            start = index
            found_Start = True
        if found_Start and amino.isalpha():
            last_amino_index = index+1

    peptide_genomique = align[0,start:last_amino_index].seq
    peptide_humain_query = align[1,start:last_amino_index].seq

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
            +str(peptide_genomique)+"\t"+str(peptide_humain_query)+"\n")
    else:
        pass
f.close()