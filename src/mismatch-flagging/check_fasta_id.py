import pandas as pd
import sys
try:
    sys.path.append('~/stage-thompson/src/')
    from Generate_Exon_Map_4 import *
except:
    print("Error import, is your stage-thompson located in home (~/) ?")
    pass

ID_file = pd.read_csv(sys.argv[1], sep="\t")
my_CDS = fasta2List(sys.argv[2])
my_genomic = fasta2List(sys.argv[3])
cds_correct = open(sys.argv[4], 'w')
genom_correct = open(sys.argv[5], 'w')

for key, value in my_CDS.items():
    my_key = key.split(" ")
    correction = ID_file.loc[ID_file["To"] == my_key[0][1:]]
    correction = correction.iloc[:, 0].to_list()
    correction = ' '.join(correction)
    cds_correct.write(my_key[0]+" "+correction+" CDS\n")
    cds_correct.write(value+"\n")

for key, value in my_genomic.items():
    my_key = key.split(" ")
    correction = ID_file.loc[ID_file["To"] == my_key[0][1:]]
    correction = correction.iloc[:, 0].to_list()
    correction = ' '.join(correction)
    genom_correct.write(my_key[0]+" "+correction+" GENOMIC\n")
    genom_correct.write(value+"\n")

cds_correct.close()
genom_correct.close()
