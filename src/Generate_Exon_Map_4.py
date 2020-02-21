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
    # Function: Merge two dictionarry together. Used to merge json response after being converted to dict.
    # Parameters:
    # 		x: (dict) first dict
    #       y: (dict) second dict
    # Return:
    # 		z: (dict) merged dict
    z = x.copy()
    z.update(y)
    return z


def importJson(filePath):
    # Function: Import json file. Json file can contains multiple json reponse with one json response to each line.
    # Parameters:
    # 		filePath: (str) path to the json file to import
    # Return:
    # 		MyJsonFull: (dict) json file merged and converted to a dictionnary
    data_dict = []
    with open(filePath) as json_data:
        for i in json_data:
            data_dict.append(json.loads(i))
    MyJsonFull = data_dict[0]
    for i in data_dict[1:]:
        MyJsonFull = merge_two_dicts(MyJsonFull, i)
    return MyJsonFull


def catExonPos(MyJsonFull, seqID, dict_Genomic):
    # Function: Extract exon position and sequence from exon informations and genomique sequence
    # Parameters:
    # 		MyJsonFull: (dict) json file merged and converted to a dictionnary with exon informations (transcript information)
    #       seqID: (str) transcript ID of the transcript of interest
    #       dict_Genomic: fasta genomic sequence converted to a dictionnary
    # Return:
    # 		exonString: (str) string with tabulated separated fields to write to a file.
    # Description: This function normalize the absolute genome position to the relative position of the genomic transcript sequence
    # depending on the strand. It also get the sequence of each exon.
    exonString = []
    try:
        if MyJsonFull[seqID]["strand"] == -1:
            ABSOLUTE_POST = MyJsonFull[seqID]["end"]
            exon_list = MyJsonFull[seqID]["Exon"]
            counter = 1
            for i in exon_list:
                stop = str(abs(i['end']-ABSOLUTE_POST))
                start = str(abs(i['start']-ABSOLUTE_POST)+1)
                res = [val for key, val in dict_Genomic.items() if seqID in key]
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
    except Exception as exc:
        pass
    return exonString


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


def writeExonmap(exon_map_File, ID_file, jsonFile, genomicFile):
    # Function: Create an exon_map file containing exon informations and sequence for each transcript of interest by calling catExonPos function.
    # Parameters:
    # 		exon_map_File: (str) path to the exon_map.tab file to write as output
    #       ID_file: (dataframe) uniprot-ensemble_transcript.tab conversion file. Column 2 correspond to ensembl ID
    #       jsonFile: (dict) json file merged and converted to a dictionnary with exon informations (transcript information)
    #       genomicFile: (dict) fasta genomic sequence converted to a dictionnary
    # Return: None. Write a tabulated exon_map file to exon_map_File path.
    with open(exon_map_File, "w") as exon_file:
        for i in ID_file.iloc[:, 1]:
            uniprot_name = ID_file.loc[ID_file["To"] == i].iloc[0, 0]
            for exon in catExonPos(jsonFile, i, genomicFile):
                exon_file.write(uniprot_name+"\t"+exon+"\n")


def writeIntronMap(intron_map_file, exon_map_File, dict_Genomic):
    # Function: Create an intron_map file containing introns informations and sequence for each transcript of interest by contrast with exon map file.
    # Parameters:
    # 		intron_map_file: (str) path & name to the intro_map file to write.
    #       exon_map_file: (str) path to the exon_map file to read to create the intron map
    # Return: None. Write a tabulated intron_map file to intron_map_file path.
    Exon_Map = pd.read_csv(exon_map_File, sep="\t", header=None)
    protein_list = list(set(Exon_Map.iloc[:, 0].to_list()))

    Intron_Map = open(intron_map_file, "w")

    for prot in protein_list:
        res = [val for key, val in dict_Genomic.items() if prot in key]
        if res == []:
            print("Error "+prot)
            continue
        subset = Exon_Map.loc[Exon_Map[0] == prot]
        intron_number = len(subset.index)-1
        for i in range(intron_number):
            intron_start = str(int(subset.iloc[i, 5]))
            intron_stop = str(int(subset.iloc[i+1, 4]))
            ensembl_ID = str(subset.iloc[i, 1])
            Intron_Map.write(prot+"\t"+ensembl_ID+"\tIntron\t" +
                             str(i+1)+"\t"+intron_start+"\t"+intron_stop+"\t"+res[0][int(intron_start):int(intron_stop)]+"\n")

    Intron_Map.close()


if __name__ == "__main__":
    ID_file = pd.read_csv(sys.argv[1], sep="\t")
    jsonFile = importJson(sys.argv[2])
    genomicFile = fasta2List(sys.argv[3])
    writeExonmap(sys.argv[4], ID_file, jsonFile, genomicFile)
    writeIntronMap(sys.argv[5], sys.argv[4], genomicFile)
