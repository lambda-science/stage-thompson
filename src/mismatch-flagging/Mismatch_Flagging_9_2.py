import pandas as pd
import sys
try:
    sys.path.append('~/stage-thompson/src/')
    from Generate_Exon_Map_4 import *
except:
    print("Error import, is your stage-thompson located in home (~/) ?")
    pass

Error_file = pd.read_csv(
    "../../data/mismatch-analysis/uniprot-error-mismatch/uniprot_new_errors_filt.txt", sep=" ", header=None)
ID_file = pd.read_csv(
    "../../data/mismatch-flagging/human_uniprot_ensembl_corrected2.tab", sep="\t")
jsonFile = importJson("../../data/mismatch-flagging/human_exon_dump.json")
genomicFile = fasta2List("../../data/mismatch-flagging/human_genomic.fasta")

my_CDS = fasta2List("../../data/mismatch-flagging/human_CDS.fasta")
exon_file = pd.read_csv(
    "../../data/mismatch-flagging/Exon_map_human.tab", sep="\t", header=None)
f = open("../../data/mismatch-flagging/mismatch_exon_pos_human.tab", "w")
f.write("Alignement\tError\tUniprotID\tPosStartError\tPosStopError\tFirstExonError\tLastExonError\n")

for index, row in Error_file.iloc[:, :].iterrows():
    fasta_name = row[0][20:-6]
    prot_name = row[2]
    error_start = row[5]
    error_stop = row[6]

    Prot_list = fasta2List("../../data/raw/uniprot-sequence/"+fasta_name)
    CDS = [val for key, val in my_CDS.items() if row[0][20:-15] in key]
    if CDS == []:
        continue
    mismtach_CDS = CDS[0][error_start*3:error_stop*3+3]
    subset = exon_file.loc[exon_file[0] == row[0][20:-15]]
    exon_number_list = subset[3].to_list()
    exon_seq_list = subset[6].to_list()

    # Interaive pop of exon list : seq. Check CDS in exon joint: TRUE = pop if False = seq important
    fini = False
    start_exon = "ERROR"
    stop_exon = "ERROR"
    while fini != True:
        exon_tuple = [(exon_number_list[i], exon_seq_list[i])
                      for i in range(len(exon_number_list))]
        for j in range(0, len(exon_tuple), 1):
            popped_exon = exon_tuple.pop(0)
            testing_condition = (''.join(mismtach_CDS) in ''.join(
                [exon_tuple[i][1] for i in range(len(exon_tuple))]))
            if testing_condition == False:
                start_exon = popped_exon[0]
                break

        exon_tuple = [(exon_number_list[i], exon_seq_list[i])
                      for i in range(len(exon_number_list))]
        for k in range(len(exon_tuple), 0, -1):
            popped_exon = exon_tuple.pop(len(exon_tuple)-1)
            testing_condition = (''.join(mismtach_CDS) in ''.join(
                [exon_tuple[i][1] for i in range(len(exon_tuple))]))
            if testing_condition == False:
                stop_exon = popped_exon[0]
                break
        f.write(row[0]+"\t"+row[1]+"\t"+row[2]+"\t"+str(row[3])+"\t" +
                str(row[4])+"\t"+str(start_exon)+"\t"+str(stop_exon)+"\n")
        fini = True
f.close()
