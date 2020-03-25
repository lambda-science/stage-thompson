import sqlite3
import pandas as pd
import os
import sys
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import math
import numpy as np


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


def openDatabseGetData(db_path):

    conn = sqlite3.connect(db_path)

    mismatch = pd.read_sql_query("SELECT * FROM mismatch", conn)
    mismatch = mismatch.astype({"exon_start_prim": "Int64", "exon_stop_prim": "Int64",
                                "exon_start_hum": "Int64", "exon_stop_hum": "Int64"})

    prot_seq_hum = pd.read_sql_query("""
    SELECT mismatch.mismatch_ID, mismatch.prot_hum, protein.sequence
    FROM mismatch
    JOIN protein ON mismatch.prot_hum = protein.prot_ID""", conn)

    prot_seq_prim = pd.read_sql_query("""
    SELECT mismatch.mismatch_ID, mismatch.prot_prim, protein.sequence
    FROM mismatch
    JOIN protein ON mismatch.prot_prim = protein.prot_ID""", conn)

    prim_exon_introns = pd.read_sql_query("""
    SELECT mismatch_ID, mismatch.prot_prim, exon_intron_map.'type', exon_intron_map.number_elem, exon_intron_map.seq
    FROM mismatch
    JOIN protein ON mismatch.prot_prim = protein.prot_ID
    JOIN transcript ON protein.transcript_ID = transcript.transcript_ID
    JOIN exon_intron_map ON transcript.transcript_ID = exon_intron_map.transcript_ID
    """, conn)

    human_exon_introns = pd.read_sql_query("""
    SELECT mismatch_ID, mismatch.prot_hum, exon_intron_map.'type', exon_intron_map.number_elem, exon_intron_map.seq
    FROM mismatch
    JOIN protein ON mismatch.prot_hum = protein.prot_ID
    JOIN transcript ON protein.transcript_ID = transcript.transcript_ID
    JOIN exon_intron_map ON transcript.transcript_ID = exon_intron_map.transcript_ID
    """, conn)

    prot_hum_iso = pd.read_sql_query(
        """SELECT * FROM protein_hum_isoform""", conn)

    return mismatch, prot_seq_hum, prot_seq_prim, prim_exon_introns, human_exon_introns, prot_hum_iso


def conserved_flag():
    # Flagging des mismatch conservés
    index_conserved = []
    for index, row in mismatch.iloc[:, :].iterrows():
        total += 1
        conserved = 0

        mySeq = fasta2List(
            "../../data/raw/uniprot-blast/results/"+row[1]+".id.fasta")
        mismtaching_seq = row[11]

        for i, j in mySeq.items():
            myAlign = pairwise2.align.localms(
                j, mismtaching_seq, 2, 0, -.5, -.5, one_alignment_only=True, score_only=True)
            if myAlign/len(mismtaching_seq) > 1.6:
                conserved += 1

        if conserved >= 4:
            index_conserved.append(index)
            counter += 1
    return index_conserved


def align_error_flag():
    # Flagging des erreurs d'alignement
    index_align_error = []
    for index, row in mismatch.iloc[:, :].iterrows():
        try:
            human_seq = prot_seq_hum.loc[prot_seq_hum["prot_hum"]
                                         == row[1]].iloc[0, 2]
        except:
            continue
        total += 1
        mismtaching_seq = row[11]
        myAlign = pairwise2.align.localms(
            human_seq, mismtaching_seq, 2, 0, -.5, -.5, one_alignment_only=True, score_only=True)
        if myAlign/len(mismtaching_seq) > 1.6:
            index_align_error.append(index)
            counter += 1
    return index_align_error


def repats_prot_flag():
    # Flagging des seq à repeats proteiq
    index_repeats = []
    for index, row in mismatch.iloc[:, :].iterrows():
        try:
            human_prot_Mismatch = prot_seq_hum.loc[prot_seq_hum["prot_hum"]
                                                   == row[1]].iloc[0, 2]
        except:
            continue

        total += 1
        human_seq = row[12]
        begin_seq = human_prot_Mismatch[:row[5]+1]
        end_seq = human_prot_Mismatch[row[6]+2:]
        human_prot_Mismatch = begin_seq + end_seq
        myAlign = pairwise2.align.localms(
            human_prot_Mismatch, human_seq, 2, 0, -.5, -.5, one_alignment_only=True, score_only=True)
        if myAlign/len(human_seq) > 1.6:
            index_repeats.append(index)
    return index_repeats


def N_genomic_flag():
    # Flagging des seq à N
    index_genom_n = []
    for index, row in mismatch.iloc[:, :].iterrows():
        my_CDS = []
        subset_exon_intron = prim_exon_introns.loc[prim_exon_introns["mismatch_ID"] == index]
        if isinstance(row[7], int) == False or isinstance(row[8], int) == False:
            continue
        try:
            for n in range(row[7], row[8]+1):
                my_exon = subset_exon_intron.loc[(subset_exon_intron["number_elem"] == n) & (
                    subset_exon_intron["type"] == "Exon")]
                my_intron = subset_exon_intron.loc[(subset_exon_intron["number_elem"] == n) & (
                    subset_exon_intron["type"] == "Intron")]
                try:
                    my_CDS.append(my_exon.iloc[0, 4])
                    my_CDS.append(my_intron.iloc[0, 4])
                except:
                    pass
            if "N" in ''.join(my_CDS):
                index_genom_n.append(index)
                counter += 1
            total += 1
        except:
            pass
    return index_genom_n


def one_vs_mult_exon_flag():
    # Flagging des 1 exon human vs multiple primate
    index_multiple_exon = []
    for index, row in mismatch.iloc[:, :].iterrows():
        if isinstance(row[7], int) == False or isinstance(row[9], int) == False:
            continue
        total += 1
        if (int(row[9])-int(row[10]) == 0) and (row[8]-row[7] >= 2):
            index_multiple_exon.append(index)
            counter += 1
    return index_multiple_exon


def non_canonical_flag():
    # Flagging des mismatch ayant un site humain de splicing non canonique
    index_non_canon = []
    for index, row in mismatch.iloc[:, :].iterrows():
        subset_intron = human_exon_introns.loc[(human_exon_introns['mismatch_ID'] == row[0]) & (
            human_exon_introns['type'] == "Intron")]
        if subset_intron.empty or row[9] == "ERROR":
            continue
        total += 1
        for i in range(int(row[9]), int(row[10])+2):
            intron = subset_intron.loc[subset_intron["number_elem"] == (i-1)]
            if intron.empty:
                continue
            intron = intron.iloc[0, :].to_list()
            # if  (intron[6][:2] == "GT" and intron[6][-2:] == "AG") \
            #    or (intron[6][:2] == "GC" and intron[6][-2:] == "AG") \
            #    or (intron[6][:2] == "AT" and intron[6][-2:] == "AC") and (len(intron[6]) > 30):
            if (intron[4][:2] == "GT" and intron[4][-2:] == "AG") and (len(intron[4]) > 30):
                pass
            else:
                index_non_canon.append(index)
                break
    return index_non_canon


def small_intron_flag():
    # Flagging mismatch with too small introns
    index_intron_small = []
    for index, row in mismatch.iloc[:, :].iterrows():
        subset_intron = prim_exon_introns.loc[(prim_exon_introns['mismatch_ID'] == row[0]) & (
            prim_exon_introns['type'] == "Intron")]
        if subset_intron.empty or row[7] == "ERROR":
            continue
        total += 1
        for i in range(int(row[7]), int(row[8])):
            intron = subset_intron.loc[subset_intron["number_elem"] == i]
            intron = intron.iloc[0, :].to_list()
            if len(intron[4]) <= 29:
                index_intron_small.append(index)
                counter += 1
                break
    return index_intron_small


def isoform_flag():
    # Isoform Flagging
    ID_file = pd.read_csv(
        "../../temp/isoform/uniprot_to_gene_human.tab", sep="\t")
    ID_transcript = pd.read_csv(
        "../../temp/isoform/gene_to_transcript.tab", sep="\t")
    ID_isoform = pd.read_csv(
        "../../temp/isoform/transcript_to_prot.tab", sep="\t")
    ID_file.rename(columns={"From": "prot", "To": "gene"}, inplace=True)
    ID_isoform.rename(columns={"From": "transcript",
                               "To": "prot_isoform"}, inplace=True)
    merged1 = ID_file.merge(ID_transcript, on="gene", how="inner")
    full_df = merged1.merge(ID_isoform, on="transcript", how="inner")

    index_isoform = []
    for index, row in mismatch.iloc[42:43, :].iterrows():
        # Récup de la list d'identifiant isoforme (bcp de traitement)
        subset = full_df.loc[full_df["prot"] == row[1]]
        subset_same = subset.loc[subset["prot_isoform"] == row[1]]
        subset_dif = subset.loc[subset["prot_isoform"] != row[1]]
        rename_prot = [row[1]+"-"+str(x+2)
                       for x in range(len(subset_same.index))]
        rename_prot = [row[1]] + rename_prot[:-1]
        subset_same.drop("prot_isoform", axis=1, inplace=True)
        subset_same["prot_isoform"] = rename_prot
        isoform_df = pd.concat([subset_same, subset_dif])
        id_isoform_hum = isoform_df["prot_isoform"].to_list()

        print(id_isoform_hum)
        # Faire les alignements pour chaque mismatch
        for i in id_isoform_hum:
            subdf = prot_hum_iso.loc[prot_hum_iso["prot_ID"] == i]
            try:
                human_seq_isoform = subdf.iloc[0, 1]
            except:
                error += 1
                continue
            if i == row[1]:
                continue
            myAlign = pairwise2.align.localms(
                human_seq_isoform, row[11], 2, 0, -.5, -.5, one_alignment_only=True, score_only=True)
            print(myAlign/len(row[11]))
            if myAlign/len(row[11]) > 1.6:
                index_isoform.append(index)
                counter += 1
                break
    return index_isoform


def dataframeFromindex(index_conserved, index_multiple_exon, index_non_canon, index_genom_n, index_intron_small, index_repeats, index_align_error, index_isoform):
    column_names = ["mismatch_ID", "conserved", "one_hum_multiple_prim", "non_canonical_hum_spl",
                    "N_in_genomic", "small_introns", "repeats_prot", "alignement_error", "human_isoform_exist"]
    df = pd.DataFrame(columns=column_names)
    df["mismatch_ID"] = mismatch["mismatch_ID"]
    df['conserved'] = np.where(df["mismatch_ID"].isin(index_conserved), 1, 0)
    df['one_hum_multiple_prim'] = np.where(
        df["mismatch_ID"].isin(index_multiple_exon), 1, 0)
    df['non_canonical_hum_spl'] = np.where(
        df["mismatch_ID"].isin(index_non_canon), 1, 0)
    df['N_in_genomic'] = np.where(df["mismatch_ID"].isin(index_genom_n), 1, 0)
    df['small_introns'] = np.where(
        df["mismatch_ID"].isin(index_intron_small), 1, 0)
    df['repeats_prot'] = np.where(df["mismatch_ID"].isin(index_repeats), 1, 0)
    df['alignement_error'] = np.where(
        df["mismatch_ID"].isin(index_align_error), 1, 0)
    df['human_isoform_exist'] = np.where(
        df["mismatch_ID"].isin(index_isoform), 1, 0)
    return df


if __name__ == "__main__":
    db_path = sys.argv[1]
    mismatch, prot_seq_hum, prot_seq_prim, prim_exon_introns, human_exon_introns, prot_hum_iso = openDatabseGetData(
        db_path)

    index_conserved = conserved_flag()
    index_multiple_exon = one_vs_mult_exon_flag()
    index_non_canon = non_canonical_flag()
    index_genom_n = N_genomic_flag()
    index_intron_small = small_intron_flag()
    index_repeats = repats_prot_flag()
    index_align_error = align_error_flag()
    index_isoform = isoform_flag()

    df = dataframeFromindex(index_conserved, index_multiple_exon, index_non_canon,
                            index_genom_n, index_intron_small, index_repeats, index_align_error, index_isoform)
    #df.to_sql(con=conn, name='mismatch_flag', index=False, if_exists="append")
