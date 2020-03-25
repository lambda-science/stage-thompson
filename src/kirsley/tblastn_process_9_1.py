# %%
# arg1: output file
# arg2: error file
# arg3: tblastn folder
import json
import pandas as pd
import sys
import os

# %%
if __name__ == "__main__":
    f = open(sys.argv[1], "w")
    f.write("Match\tPrimate\tHuman\tSeq Primate\tSeq Human\tStart Genome\t" +
            "Stop Genome\tE-value\tSimilarity\tLength\tFrame\n")

    Error_file = pd.read_csv(sys.argv[2], sep=" ", header=None)

    for index, row in Error_file.iloc[:, :].iterrows():
        fasta_name = row[0][58:-6]
        prot_name = row[2].split("_")
        prot_name = prot_name[0]
        human_start = row[5]
        human_stop = row[6]

        blast_file = sys.argv[3]+"blast_out/" + \
            str(index)+"_" + str(row[0][58:-30])+"_"+str(prot_name)+".blast"
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
