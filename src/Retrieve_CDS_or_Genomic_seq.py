# Script créer pour récupérer de façon parallèle les ficier des CDS des identifiants uniprot à mismatch.
# (sous-produit 3-Retrieve_Genomic....)
import pandas as pd
import grequests
import json

ID_file = pd.read_csv("transcript_ensembl.tab", sep="\t")

url = "https://rest.ensembl.org/sequence/id"
headers = {"Content-Type": "application/json", "Accept": "application/json"}

params = []
for i in range(0, 9960, 50):
    try:
        params.append({"ids": ID_file.iloc[i:i+50, 1].tolist(), "type": "cds"})
    except:
        params.append({"ids": ID_file.iloc[i:, 1].tolist(), "type": "cds"})

rs = [grequests.post(url, headers=headers, data=json.dumps(i)) for i in params]
all_response = grequests.map(rs, size=20)
print(all_response)

f = open("cds_new.fa", "w")
indice = 0
for response in all_response:
    j = 0
    for entry in response.json():
        f.write(">"+entry['query']+" " +
                str(ID_file.iloc[indice*50+j, 0]) + " CDS\n")
        f.write(entry['seq']+"\n")
        j = j+1
    indice = indice + 1
f.close()
