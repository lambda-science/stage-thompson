# Script créer pour récupérer de façon parallèle les fichier JSON contenant l'information d'exon-map
# (sous-produit 3-Retrieve_Genomic....)
import grequests
import requests
import sys
import pandas as pd
import json

ID_file = pd.read_csv("transcript_ensembl.tab", sep="\t")
url = "https://rest.ensembl.org/lookup/id/"
headers = {"Content-Type": "application/json", "Accept": "application/json"}

params = []
for i in range(0, 9960, 1000):
    try:
        params.append(
            {"ids": ID_file.iloc[i:i+1000, 1].tolist(), "expand": "1"})
    except:
        params.append({"ids": ID_file.iloc[i:, 1].tolist(), "expand": "1"})

rs = [grequests.post(url, headers=headers, data=json.dumps(i)) for i in params]
all_response = grequests.map(rs)

print(all_response)
f = open("../raw/uniprot-exon-map/lookup_newline.json", "a")
for response in all_response:
    try:
        f.write(json.dumps(response.json())+"\n")
    except:
        pass
f.close()
