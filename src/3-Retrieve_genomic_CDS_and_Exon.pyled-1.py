# # Script permettant de récupérer les séquences genomiques/CDS des transcript d'intéret (sequence/id) ainsi que les positiosn des exons (lookup/id)
#
# Input:
# arg1 : fichier transcript_ensembl.tab tableau de correspondance Uniprot AC/ID -> Ensembl Transcript téléchargé sur https://www.uniprot.org/uploadlists/
# arg2: fichier ou écrire les CDS
# arg3: ficheir ou écrire les sequences GENOMIQUES
# arg4: fichier ou écrire l'exon MAP

# Output:
#     genomic_new.fa / cds_new.fa fichier contenant toute les séquences fasta (génomique ou CDS) pour nos protéines d'intérêt
#     lookup_newline.json fichier JSON contenant les coordonées exonique de toute les séquences à traiter pour générer les exon map.

# Importation des lib et data

import pandas as pd
import requests
import sys
import grequests
import json


def makeAsyncEnsemblSeqRequest(ID_file, type_request):
    # Création des requêtes parallèles sur les serveur de ensembl pour récupérer les séquences génomique ou CDS
    # Penser à indiquer le bon type: genomic ou cds
    # Peut etre assez long (1h - 1H30), changer le paramètre size pour faire plus de requête simultanée
    # Check en printant all_response si toute les réponses sont 200, dans le cas contraire certaines données sont perdues / sont à refaire
    url = "https://rest.ensembl.org/sequence/id"
    headers = {"Content-Type": "application/json",
               "Accept": "application/json"}

    params = []
    for i in range(0, len(ID_file.index), 50):
        try:
            params.append(
                {"ids": ID_file.iloc[i:i+50, 1].tolist(), "type": type_request})
        except:
            params.append(
                {"ids": ID_file.iloc[i:, 1].tolist(), "type": type_request})

    rs = [grequests.post(url, headers=headers, data=json.dumps(i))
          for i in params]
    all_response = grequests.map(rs)
    return all_response


def writeAsyncEnsemblResponse(all_response, ID_file, file_path, type_request):
    # Ecriture des résultats des réponses dans un fichier fasta. Penser à indiquer le bon type CDS ou GENOMIC
    f = open(file_path, "w")
    indice = 0
    for response in all_response:
        j = 0
        try:
            for entry in response.json():
                f.write(
                    ">"+entry['query']+" "+str(ID_file.iloc[indice*50+j, 0]) + type_request + "\n")
                f.write(entry['seq']+"\n")
                j = j+1
        except:
            pass
        indice = indice + 1
    f.close()


def makeAsyncEnsemblExonmapRequest(ID_file):
    # Création des requêtes parallèles sur les serveur de ensembl pour récupérer les séquences informations de localisation des exons
    # Peut etre assez long (1h)
    # Check en printant all_response si toute les réponses sont 200, dans le cas contraire certaines données sont perdues / sont à refaire
    url = "https://rest.ensembl.org/lookup/id/"
    headers = {"Content-Type": "application/json",
               "Accept": "application/json"}

    params = []
    for i in range(0, len(ID_file.index), 1000):
        try:
            params.append(
                {"ids": ID_file.iloc[i:i+1000, 1].tolist(), "expand": "1"})
        except:
            params.append({"ids": ID_file.iloc[i:, 1].tolist(), "expand": "1"})

    rs = [grequests.post(url, headers=headers, data=json.dumps(i))
          for i in params]
    all_response = grequests.map(rs)
    return all_response


def writeAsyncEnsemblExonMapResposne(all_response, file_write):
    # Ecriture des résultats des réponses dans un fichier json.
    f = open(file_write, "a")
    for response in all_response:
        try:
            f.write(json.dumps(response.json())+"\n")
        except:
            pass
    f.close()


if __name__ == "__main__":
    #ID_file = pd.read_csv("../raw/uniprot-exon-map/transcript_ensembl.tab", sep = "\t")
    ID_file = pd.read_csv(sys.argv[1], sep="\t")

    my_response = makeAsyncEnsemblSeqRequest(ID_file, "cds")
    print(my_response)
    writeAsyncEnsemblResponse(my_response, ID_file, sys.argv[2], "CDS")

    my_response2 = makeAsyncEnsemblSeqRequest(ID_file, "genomic")
    print(my_response2)
    writeAsyncEnsemblResponse(my_response2, ID_file, sys.argv[3], "GENOMIC")

    my_response3 = makeAsyncEnsemblExonmapRequest(ID_file)
    print(my_response3)
    writeAsyncEnsemblExonMapResposne(my_response3, sys.argv[4])
