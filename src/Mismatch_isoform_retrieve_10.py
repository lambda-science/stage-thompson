import sqlite3
import pandas as pd
import os
import sys
import grequests
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import math
import numpy as np
import requests
import json


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


def prot_to_gene():
    # Flagging isoform
    # Prot to gene
    f = open("../temp/isoform/uniprot_to_gene_human.tab", "w")
    uniq_human_prot = list(set(mismatch["prot_hum"].to_list()))
    URL = "https://www.uniprot.org/uploadlists/"
    params = {
        "from": "ACC+ID",
        "to": "ENSEMBL_ID",
        "format": "tab",
        "query": ' '.join(uniq_human_prot)
    }

    r = requests.post(URL, data=params)
    f.write(r.text)
    f.close()


def gene_to_transcript_json():
    # Flagging isoform
    # Gene to Transcript
    error = True
    ID_file = pd.read_csv(
        "../temp/isoform/uniprot_to_gene_human.tab", sep="\t")
    url = "https://rest.ensembl.org/lookup/id/"
    headers = {"Content-Type": "application/json",
               "Accept": "application/json"}

    params = []
    for i in range(0, len(ID_file.index), 300):
        try:
            params.append(
                {"ids": ID_file.iloc[i:i+300, 1].tolist(), "expand": "1"})
        except:
            params.append({"ids": ID_file.iloc[i:, 1].tolist(), "expand": "1"})

    rs = [grequests.post(url, headers=headers, data=json.dumps(i))
          for i in params]
    all_response = grequests.map(rs, size=3)

    while error == True:
        error = False
        for index, response in enumerate(all_response):
            if response is None:
                error = True
                r = requests.post(url, headers=headers,
                                  data=json.dumps(params[index]))
                all_response[index] = r
                continue

            elif not response.ok:
                error = True
                r = requests.post(url, headers=headers,
                                  data=json.dumps(params[index]))
                all_response[index] = r

    f = open("../temp/isoform/json_dump_transcript.json", "a")
    for response in all_response:
        try:
            f.write(json.dumps(response.json())+"\n")
        except:
            pass
    f.close()


def write_gene_to_transcript():
    transcript = importJson("../temp/isoform/json_dump_transcript.json")
    f = open("../temp/isoform/gene_to_transcript.tab", "w")
    f.write("gene\ttranscript\n")
    for i in transcript:
        for j in range(len(transcript[i]["Transcript"])):
            f.write(i+"\t"+transcript[i]["Transcript"][j]["id"]+"\n")
    f.close()


def transcript_to_prot():
    # Transcript to prot
    # Flagging isoform

    f = open("../temp/isoform/transcript_to_prot.tab", "w")

    ID_transcript = pd.read_csv(
        "../temp/isoform/gene_to_transcript.tab", sep="\t")
    uniq_human_trans = list(set(ID_transcript["transcript"].to_list()))

    URL = "https://www.uniprot.org/uploadlists/"
    params = {
        "from": "ENSEMBL_TRS_ID",
        "to": "ACC",
        "format": "tab",
        "query": ' '.join(uniq_human_trans)
    }

    r = requests.post(URL, data=params)
    f.write(r.text)
    f.close()


def write_protein_iso_db():
    # Dataframe pour Protein
    uniprot_ID = []
    uniprot_Seq = []
    Prot_list = fasta2List("../temp/isoform/all_isoform.id.fasta")
    for key, val in Prot_list.items():
        myKey = key[1:].split(" ")
        uniprot_ID.append(myKey[0])
        uniprot_Seq.append(val)
    dict_uniprot = {"uniprot_ID": uniprot_ID, "uniprot_Seq": uniprot_Seq}
    df_prot = pd.DataFrame(dict_uniprot)
    df_prot["organism"] = 9606
    df_prot = df_prot.rename(
        columns={"uniprot_ID": "prot_ID", "uniprot_Seq": "sequence"})
    df_prot.to_sql(con=conn, name='protein_hum_isoform',
                   index=False, if_exists="append")


if __name__ == "__main__":
    db_path = sys.argv[1]
    conn = sqlite3.connect(db_path)

    mismatch = pd.read_sql_query("SELECT * FROM mismatch", conn)
    mismatch = mismatch.astype({"exon_start_prim": "Int64", "exon_stop_prim": "Int64",
                                "exon_start_hum": "Int64", "exon_stop_hum": "Int64"})

    prot_to_gene()
    gene_to_transcript_json()
    write_gene_to_transcript()
    transcript_to_prot()
    # write_protein_iso_db()
