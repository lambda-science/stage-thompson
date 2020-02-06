# Script permettant de générer la carte exonique (exon map) des gènes d'intérêt
# Input:
# arg1: ID_file (correspondance uniprot-ensembl tabulée transcript_ensembl.tab)
# arg2: jsonFile (json généré par script 3 des positions exoniques)
# arg3: my_Genomic (fichier contenant les séquences génomiques au format fasta)
# arg4: exon_map_file (fichier d'écriture de l'exon map)
# arg5: intron_map_file (fichier d'écrite de l'intron map

# Output: exon_map.tab / intron_map fichier texte tabulation contenant la position de chaque exon et la séquence pour chaque transcript d'intérêt

import json
import pandas as pd
import sys


def merge_two_dicts(x, y):
    # Pour merger les json entre eux
    z = x.copy()
    z.update(y)
    return z


def importJson(filePath):
    # Pour importer les json
    data_dict = []
    with open(filePath) as json_data:
        for i in json_data:
            data_dict.append(json.loads(i))
    MyJsonFull = data_dict[0]
    for i in data_dict[1:]:
        MyJsonFull = merge_two_dicts(MyJsonFull, i)
    return MyJsonFull


def catExonPos(MyJsonFull, seqID, dict_Genomic):
    # Renvoie une string contenant la ligne à écrire dans le fichier (position exon, séquence...)
    exonString = []
    if MyJsonFull[seqID]["strand"] == -1:
        ABSOLUTE_POST = MyJsonFull[seqID]["end"]
        exon_list = MyJsonFull[seqID]["Exon"]
        counter = 1
        for i in exon_list:
            stop = str(abs(i['end']-ABSOLUTE_POST))
            start = str(abs(i['start']-ABSOLUTE_POST)+1)
            res = [val for key, val in my_Genomic.items() if seqID in key]
            exonString.append(seqID + "\tExon\t"+str(counter)+"\t" +
                              stop+"\t"+start+"\t"+res[0][int(stop):int(start)])
            counter += 1

    if MyJsonFull[seqID]["strand"] == 1:
        ABSOLUTE_POST = MyJsonFull[seqID]["start"]
        exon_list = MyJsonFull[seqID]["Exon"]
        counter = 1
        for i in exon_list:
            start = str(abs(i['start']-ABSOLUTE_POST))
            stop = str(abs(i['end']-ABSOLUTE_POST)+1)
            res = [val for key, val in dict_Genomic.items() if seqID in key]
            exonString.append(seqID + "\tExon\t"+str(counter)+"\t" +
                              start+"\t"+stop+"\t"+res[0][int(start):int(stop)])
            counter += 1
    return exonString


def fasta2List(pathFasta):
    # Importe un fichier fasta en dictionnaire titre:séquence
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


def writeExonmap(exon_map_File, ID_file, jsonFile, genomicFile):
    with open("../raw/uniprot-exon-map/Exon_map.tab", "w") as exon_file:
        for i in ID_file.iloc[:, 1]:
            uniprot_name = ID_file.loc[ID_file["To"] == i].iloc[0, 0]
            for exon in catExonPos(jsonFile, i, genomicFile):
                exon_file.write(uniprot_name+"\t"+exon+"\n")


def writeIntronMap(intron_map_file, exon_map_File):
    # On fait maintenant l'intron map pour trouver les petits introns
    Exon_Map = pd.read_csv(exon_map_File, sep="\t", header=None)
    protein_list = list(set(Exon_Map.iloc[:, 0].to_list()))

    Intron_Map = open(intron_map_file, "w")

    for prot in protein_list:
        subset = Exon_Map.loc[Exon_Map[0] == prot]
        intron_number = len(subset.index)-1
        for i in range(intron_number):
            intron_start = str(int(subset.iloc[i, 5]))
            intron_stop = str(int(subset.iloc[i+1, 4]))
            ensembl_ID = str(subset.iloc[i, 1])
            Intron_Map.write(prot+"\t"+ensembl_ID+"\tIntron\t" +
                             str(i+1)+"\t"+intron_start+"\t"+intron_stop+"\n")

    Intron_Map.close()


if __name__ == "__main__":
    ID_file = pd.read_csv(sys.argv[1], sep="\t")
    jsonFile = importJson(sys.argv[2])
    my_Genomic = fasta2List(sys.argv[3])
    writeExonmap(sys.argv[4], ID_file, jsonFile, my_Genomic)
    writeIntronMap(sys.argv[5], sys.argv[4])
