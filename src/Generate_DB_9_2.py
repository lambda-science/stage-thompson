import numpy as np
import sys
import os
import pandas as pd
import sqlite3
conn = sqlite3.connect('../../mismatch_db.db')


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


# Table mismatch
error_list = pd.read_csv(
    "../../data/mismatch-analysis/uniprot_errors_type3.txt", sep=" ", header=None)
human_exon_pos = pd.read_csv(
    "../../data/mismatch-flagging/human_mismatch_exon_pos.tab", sep="\t")
primate_exon_pos = pd.read_csv(
    "../../data/mismatch-analysis/mismatch_exon_pos.tab", sep="\t")

# Script kirsley pr seq
mismatch_prim = []
mismatch_hum = []

for index, row in error_list.iloc[:, :].iterrows():
    fasta_name = row[0][58:-6]
    Prot_list = fasta2List("../../data/raw/uniprot-sequence/"+fasta_name)
    prot_prim = [val for key, val in Prot_list.items() if row[2] in key]
    prot_hum = [val for key, val in Prot_list.items() if row[0]
                [58:-15] in key]
    mismatch_prim.append(prot_prim[0][row[3]:row[4]+1])
    mismatch_hum.append(prot_hum[0][row[5]:row[6]+1])
error_list.insert(7, "seq_prim", mismatch_prim, allow_duplicates=True)
error_list.insert(8, "seq_hum", mismatch_hum, allow_duplicates=True)

# Rename des col
primate_exon_pos = primate_exon_pos.rename({"Alignement": "prot_hum", "UniprotID": "prot_prim", "PosStartError": "pos_start_prim",
                                            "PosStopError": "pos_stop_prim", "FirstExonError": "exon_start_prim", "LastExonError": "exon_stop_prim"}, axis=1)
human_exon_pos.drop(["PosStartError_H", "PosStopError_H"],
                    axis=1, inplace=True)
human_exon_pos = human_exon_pos.rename({"Alignement": "prot_hum", "UniprotID": "prot_prim", "PosStartError_P": "pos_start_prim",
                                        "PosStopError_P": "pos_stop_prim", "FirstExonError": "exon_start_hum", "LastExonError": "exon_stop_hum"}, axis=1)
error_list = error_list.rename({0: "prot_hum", 1: "Error", 2: "prot_prim", 3: "pos_start_prim",
                                4: "pos_stop_prim", 5: "pos_start_hum", 6: "pos_stop_hum"}, axis=1)
# Remplacement des erreur par NAN
primate_exon_pos.replace("ERROR", np.nan, inplace=True)
human_exon_pos.replace("ERROR", np.nan, inplace=True)
# Merging des trois tableaux
exon_pos = human_exon_pos.merge(primate_exon_pos, how="outer", on=[
                                "prot_hum", "Error", "prot_prim", "pos_start_prim", "pos_stop_prim"])
error_list_db = error_list.merge(exon_pos, how='outer', on=[
                                 "prot_hum", "Error", "prot_prim", "pos_start_prim", "pos_stop_prim"])
# Retyping (bug merge)
error_list_db = error_list_db.astype(
    {"exon_start_prim": "float", "exon_stop_prim": "float", "exon_start_hum": "float", "exon_stop_hum": "float"})
error_list_db = error_list_db.astype(
    {"exon_start_prim": "Int64", "exon_stop_prim": "Int64", "exon_start_hum": "Int64", "exon_stop_hum": "Int64"})
# Reordering and substring
error_list_db["mismatch_ID"] = error_list_db.index
error_list_db = error_list_db[["mismatch_ID", "prot_hum", "prot_prim", "pos_start_prim", "pos_stop_prim", "pos_start_hum",
                               "pos_stop_hum", "exon_start_prim", "exon_stop_prim", "exon_start_hum", "exon_stop_hum", "seq_prim", "seq_hum"]]
error_list_db['prot_hum'] = error_list_db['prot_hum'].str[58:-15]
error_list_db.to_sql(con=conn, name='mismatch',
                     index=False, if_exists="append")

# Table Protein
primate_prot = fasta2List("../../data/raw/uniprot-sequence/all_sequence.fasta")
human_prot = fasta2List("../../data/raw/uniprot-sequence/all_sequence.fasta")
primate_ensembl = pd.read_csv(
    "../../data/mismatch-analysis/transcript_ensembl_corrected2.tab", sep="\t")
human_ensembl = pd.read_csv(
    "../../data/mismatch-flagging/human_transcript_ensembl_corrected2.tab", sep="\t")

primate_prot.update(human_prot)
ensembl_uniprot = pd.concat([human_ensembl, primate_ensembl])
ensembl_uniprot.rename(
    columns={'From': "prot_ID", 'To': "transcript_ID"}, inplace=True)
prot_name = []
orga = []
seq = []

for key, val in primate_prot.items():
    my_elem = key.split(" ")
    prot_name.append(my_elem[0][1:])
    orga.append([x for x in my_elem if "OX=" in x][0][3:])
    seq.append(val)
my_dict = {'prot_ID': prot_name, 'sequence': seq, 'organism': orga}
df = pd.DataFrame(my_dict)
prot_table = df.merge(ensembl_uniprot, on=["prot_ID"], how='left')
prot_table.to_sql(con=conn, name='protein', index=False, if_exists="append")

# Table Transcript
primate_ensembl = pd.read_csv(
    "../../data/mismatch-analysis/transcript_ensembl_corrected2.tab", sep="\t")
human_ensembl = pd.read_csv(
    "../../data/mismatch-flagging/human_transcript_ensembl_corrected2.tab", sep="\t")
primate_cds = fasta2List("../../data/mismatch-analysis/CDS_all_filt.fasta")
human_cds = fasta2List(
    "../../data/mismatch-flagging/human_CDS_all_corrected_filt.fasta")
primate_genom = fasta2List(
    "../../data/mismatch-analysis/genomic_all_filt.fasta")
human_genom = fasta2List(
    "../../data/mismatch-flagging/human_genomic_all_corrected_filt.fasta")
human_cds.update(primate_cds)
human_genom.update(primate_genom)
ensembl_uniprot = pd.concat([human_ensembl, primate_ensembl])

# Manipulation, merging, écriture
df_cds = pd.DataFrame({'transcript_ID': list(
    human_cds.keys()), 'sequence_CDS': list(human_cds.values())})
df_cds['transcript_ID'] = df_cds['transcript_ID'].str.split(" ")
df_cds['transcript_ID'] = df_cds['transcript_ID'].str[0]
df_cds['transcript_ID'] = df_cds['transcript_ID'].str[1:]
df_genom = pd.DataFrame({'transcript_ID': list(
    human_genom.keys()), 'sequence_genomic': list(human_genom.values())})
df_genom['transcript_ID'] = df_genom['transcript_ID'].str.split(" ")
df_genom['transcript_ID'] = df_genom['transcript_ID'].str[0]
df_genom['transcript_ID'] = df_genom['transcript_ID'].str[1:]
df_cds = df_cds.merge(df_genom, on='transcript_ID', how='outer')
df_cds.to_sql(con=conn, name='transcript', index=False, if_exists="append")

# Table exon_intron_map
# stage-thompson/data$ cat mismatch-analysis2/Exon_map.tab mismatch-analysis2/Intron_map.tab mismatch-flagging/human_Exon_map.tab mismatch-flagging/human_Intron_map.tab > full_exon_intron_map.tab
full_map = pd.read_csv(
    "../../data/full_exon_intron_map.tab", sep="\t", header=None)
full_map.drop(0, axis=1, inplace=True)
full_map = full_map.rename({1: "transcript_ID", 2: "type", 3: "number_elem",
                            4: "pos_start_genom", 5: "pos_end_genom", 6: "seq"}, axis=1)
full_map["transcript_index"] = full_map.index
full_map = full_map[["transcript_index", "transcript_ID", "type",
                     "number_elem", "pos_start_genom", "pos_end_genom", "seq"]]
full_map.to_sql(con=conn, name='exon_intron_map',
                index=False, if_exists="append")

# Table tblastn
tblastn = pd.read_csv(
    "../../data/mismatch-analysis/tblastn/match.out", sep="\t")
tblastn.drop(["Match", "Similarity", "Length", "Frame"], axis=1, inplace=True)
tblastn = tblastn.rename({"Primate": "prot_ID_prim", "Human": "prot_ID_hum", "Seq Primate": "seq_in_prim",
                          "Seq Human": "peptide_hum", "Start Genome": "start_genom", "Stop Genome": "stop_genom", "E-value": "e_value"}, axis=1)
tblastn.to_sql(con=conn, name='tblastn_match', index=False, if_exists="append")
