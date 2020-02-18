# IMPORTANT: you need to  meyer@milex:~/stage-thompson/data/raw/uniprot-sequence$ cat *.fasta > all_sequence.fasta
# in order to be able to do that
import pandas as pd
import sys
try:
    sys.path.append('~/stage-thompson/src/')
    from Generate_Exon_Map_4 import *
except:
    pass

ID_file = pd.read_csv(
    "../data/mismatch-analysis/uniprot-exon-map/transcript_ensembl.tab", sep="\t")
my_CDS = fasta2List(
    "../data/mismatch-analysis/uniprot-exon-map/cds_new.fa")
Prot_list = fasta2List(
    "../data/raw/uniprot-sequence/all_sequence.fasta")

f = open("../data/mismatch-analysis/uniprot-exon-map/transcript_ensembl_corrected.tab", "w")
f.write("From\tTo\n")

for index, row in ID_file.iloc[:, :].iterrows():

    CDS = [val for key, val in my_CDS.items() if row[1] in key]
    prot = [val for key, val in Prot_list.items() if row[0] in key]

    if CDS == [] or prot == []:
        continue
    elif len(CDS[0])-3 == len(prot[0])*3:
        f.write(row[0] + "\t" + row[1]+"\n")

f.close()
