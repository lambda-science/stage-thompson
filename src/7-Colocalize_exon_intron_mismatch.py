# Script permettant de co-localiser les mismatch au niveau introns / exons
# INPUT:
# arg1: error_file: fichier list d'erreur d'alignement (script julie)
# arg2: my_CDS: CDS des séquence présentant une erreur d'intérêt
# arg3: exon_file: exon map des séquences présentant une erreur d'intérêt
# arg4: mismatch_exon_pos_file: fichier (output 1ere fonction) de localisation du mismatch au niveau des exons
# arg5: intron_file: intron map des séquences présentant une erreur d'intérêt
# arg6: exon_out: chemin fichier de sortie de la séquences des introns (mismatch localisé)
# arg7: intron_out: chemin fichier de sortie de la séquences des exons (mismatch localisé)
# OUTPUT:
# mismatch_exon_pos_file, exon_out, intron_out

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

def localizeMismatchOnExonMap(mismatch_exon_pos_file, Error_file, my_CDS, exon_file):
    # Génération du fichier qui map chaque mismatch à un exon de début et de fin de mismatch
    f = open(mismatch_exon_pos_file, "w")
    f.write("Alignement\tError\tUniprotID\tPosStartError\tPosStopError\tFirstExonError\tLastExonError\n")
    for index, row in Error_file.iloc[:,:].iterrows():
        prot_name = row[2]
        error_start = row[3]
        error_stop = row[4]

        CDS = [val for key, val in my_CDS.items() if prot_name in key]
        if CDS == []:
            continue
        mismtach_CDS = CDS[0][error_start*3:error_stop*3+3]
        subset = exon_file.loc[exon_file[0]==prot_name]
        exon_number_list = subset[3].to_list()
        exon_seq_list = subset[6].to_list()
        
        # On enlève des exons aux extrémités de façon itérative
        # Si le mismatch ne se trouve plus dans la liste d'exon
        # Alors cela veut dire que l'exon supprimé contenait le mismatch
        fini = False
        while fini != True:
            exon_tuple = [ (exon_number_list[i], exon_seq_list[i]) for i in range(len(exon_number_list))]
            for j in range(0,len(exon_tuple), 1):
                popped_exon = exon_tuple.pop(0)
                testing_condition = (''.join(mismtach_CDS) in ''.join( [ exon_tuple[i][1] for i in range(len(exon_tuple)) ]))
                if testing_condition == False:
                    start_exon = popped_exon[0]
                    break

            exon_tuple = [ (exon_number_list[i], exon_seq_list[i]) for i in range(len(exon_number_list))]
            for k in range(len(exon_tuple), 0, -1):
                popped_exon = exon_tuple.pop(len(exon_tuple)-1)
                testing_condition = (''.join(mismtach_CDS) in ''.join( [ exon_tuple[i][1] for i in range(len(exon_tuple)) ]))
                if testing_condition == False:
                    stop_exon = popped_exon[0]
                    break
            f.write(row[0]+"\t"+row[1]+"\t"+row[2]+"\t"+str(row[3])+"\t"+str(row[4])+"\t"+str(start_exon)+"\t"+str(stop_exon)+"\n")
            fini = True
    f.close()

def getExonIntronMismatchSeq(mismatch_exon_pos_file, exon_file, intron_file, exon_out, intron_out):
    # Importation de tout les introns / exons présent dans les mismatch pour regarder la distribution des tailles
    mismatch_pos = pd.read_csv(mismatch_exon_pos_file, sep="\t")
    # Récupération des séquences des exons introns localisés dans la zone du mismatch
    f = open(exon_out, "w", newline='\n')
    f2 = open(intron_out, "w", newline='\n')

    for index, row in mismatch_pos.iloc[:,:].iterrows():
        subset_exon = exon_file.loc[exon_file[0]==row[2]]
        subset_intron = intron_file.loc[intron_file[0]==row[2]]
        for i in range(row[5], row[6]+1):
            row_to_list = subset_exon.loc[subset_exon[3] == i].values.tolist()
            row_to_list = row_to_list[0]
            my_Str = '\t'.join(map(str, row_to_list))
            f.write(my_Str+"\n")
        for i in range(row[5], row[6]):
            row_to_list = subset_intron.loc[subset_intron[3] == i].values.tolist()
            row_to_list = row_to_list[0]
            my_Str = '\t'.join(map(str, row_to_list))
            f2.write(my_Str+"\n")
            pass
    f.close()
    f2.close()

# Importation de toutes les données utilisées
if __name__=="__main__":
    Error_file = pd.read_csv(sys.argv[1], sep=" ", header=None)
    my_CDS = fasta2List(sys.argv[2])
    exon_file = pd.read_csv(sys.argv[3], sep="\t", header=None)
    mismatch_exon_pos_file = sys.argv[4]
    
    localizeMismatchOnExonMap(mismatch_exon_pos_file, Error_file, my_CDS, exon_file)

    intron_file = pd.read_csv(sys.argv[5], sep="\t", header=None)
    exon_out = sys.argv[6]
    intron_out = sys.argv[7]
    getExonIntronMismatchSeq(mismatch_exon_pos_file, exon_file, intron_file, exon_out, intron_out)