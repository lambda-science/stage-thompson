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


def makeAsyncEnsemblExonmapRequest(ID_file):
    # Function: Create async request on Ensembl API to get exon informations about a list of transcripts.
    # Parameters:
    # 		ID_file: (dataframe) uniprot-ensemble_transcript.tab conversion file. Column 2 correspond to ensembl ID
    # Return:
    # 		all_response: (list) list of response to ensemble API request (json format)
    # Description: This function can have a long running time (multiple hours) you can change the size=10 parameters to increase
    # the number of simultaneous request. You can print all_reposne to check if all response are 200 meaning they all worked.
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
    all_response = grequests.map(rs, size=1)
    return all_response


def writeAsyncEnsemblExonMapResposne(all_response, file_write):
    # Function: Write list of response json results to a single json file.
    # Parameters:
    # 		all_response: (list) list of response to ensemble API request (json format)
    #       file_write: (str) path to the json file to write
    # Return: None - Write a file at file_write
    f = open(file_write, "a")
    for response in all_response:
        try:
            f.write(json.dumps(response.json())+"\n")
        except:
            pass
    f.close()


if __name__ == "__main__":
    ID_file = pd.read_csv(sys.argv[1], sep="\t")

    my_response3 = makeAsyncEnsemblExonmapRequest(ID_file)
    print(my_response3)
    writeAsyncEnsemblExonMapResposne(my_response3, sys.argv[2])
