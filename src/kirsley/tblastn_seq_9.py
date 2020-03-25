# %%
import json
import pandas as pd
import sys
import os
# %%


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


def writeQuerySubject(Error_file, out_folder, my_Genomic, raw_folder):
    os.mkdir(out_folder+"/query_subject")
    os.mkdir(out_folder+"/blast_out")

    fff = open(out_folder+"/all_couple.txt", "w")
    for index, row in Error_file.iloc[:, :].iterrows():
        fasta_name = row[0][58:-6]
        prot_name = row[2].split("_")
        prot_name = prot_name[0]

        human_start = row[5]
        human_stop = row[6]
        Prot_list = fasta2List(raw_folder+"/"+fasta_name)

        prot_HumanRef = [val for key, val in Prot_list.items()
                         if row[0][58:-15] in key]
        genomic_Seq = [val for key, val in my_Genomic.items()
                       if prot_name in key]
        if genomic_Seq == []:
            pass
        elif prot_HumanRef == []:
            pass

        else:
            peptide_Ref = (row[0][58:-15], prot_HumanRef[0]
                           [human_start:human_stop])
            f = open(out_folder+"/query_subject/"+str(index)+"_" +
                     str(row[0][58:-15])+"_"+str(prot_name)+".subject", "w")
            ff = open(out_folder+"/query_subject/"+str(index)+"_" +
                      str(row[0][58:-15])+"_"+str(prot_name)+".query", "w")

            f.write(">"+prot_name+"\n"+str(genomic_Seq[0]))
            ff.write(">"+row[0][58:-15]+"\n"+str(peptide_Ref[1]))
            fff.write(str(index)+"_" +
                      str(row[0][58:-15])+"_"+str(prot_name)+"\n")

            f.close()
            ff.close()
    fff.close()


if __name__ == "__main__":
    Error_file = pd.read_csv(sys.argv[1], sep=" ", header=None)
    my_Genomic = fasta2List(sys.argv[2])
    raw_folder = sys.argv[3]
    out_folder = sys.argv[4]
    writeQuerySubject(Error_file, out_folder, my_Genomic, raw_folder)
