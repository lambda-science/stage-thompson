# %%
import json
import pandas as pd
import sys
import os

# %%


def processtBlastn(Error_file, blast_folder):
    for index, row in Error_file.iloc[:, :].iterrows():
        fasta_name = row[0][20:-6]
        prot_name = row[2]
        human_start = row[5]
        human_stop = row[6]

        blast_file = blast_folder+"/" + \
            str(index)+"_" + str(row[0][20:-15])+"_"+str(prot_name)+".blast"
        df = pd.read_csv(blast_file, sep="\t", header=None)


# %%
if __name__ == "__main__":
    f = open("data/mismatch-analysis/tblastn/match.out", "w")
    f.write("Match\tPrimate\tHuman\tSeq Primate\tSeq Human\tStart Genome\t" +
            "Stop Genome\tE-value\tSimilarity\tLength\tFrame\n")
    Error_file = pd.read_csv(
        "data/mismatch-analysis/uniprot-error-mismatch/uniprot_new_errors_" +
        "filt.txt", sep=" ", header=None)
    for index, row in Error_file.iloc[:, :].iterrows():
        fasta_name = row[0][20:-6]
        prot_name = row[2]
        human_start = row[5]
        human_stop = row[6]

        blast_file = "data/mismatch-analysis/tblastn/blast_out/" + \
            str(index)+"_" + str(row[0][20:-15])+"_"+str(prot_name)+".blast"
        try:
            df = pd.read_csv(blast_file, sep="\t", header=None)
        except:
            pass
        if 'df' in locals():
            match = df.iloc[0, :].tolist()
            if match[7] / match[8] >= 0.90 and match[6] < 0.001 \
                    and (match[9] not in [-1, -2, -3]):
                line = "Match\t"+"\t".join(str(x) for x in match)+"\n"
                f.write(line)
            del(df)
    f.close()

# %%
