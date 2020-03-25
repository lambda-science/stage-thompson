# arg1: error_file
# arg2: error_file out
import pandas as pd
try:
    import sys
    sys.path.append('~/stage-thompson/src/')
    from Generate_Exon_Map_4 import fasta2List
except:
    pass

if __name__ == "__main__":
    Error_file = pd.read_csv(sys.argv[1], sep=" ", header=None)
    f = open(sys.argv[2], 'w')
    f.write("Alignement_File Error_Type UniprotID StartPrimate StopPrimate StartHuman StopHuman SeqPrimate SeqHuman\n")
    for index, row in Error_file.iloc[:, :].iterrows():
        fasta_name = row[0][58:-6]
        Prot_list = fasta2List("../data/raw/kirsley/"+fasta_name)
        prot_prim = [val for key, val in Prot_list.items() if row[2] in key]
        prot_hum = [val for key, val in Prot_list.items() if row[0]
                    [58:-30] in key]
        mismatch_prim = prot_prim[0][row[3]:row[4]+1]
        mismatch_hum = prot_hum[0][row[5]:row[6]+1]

        f.write(' '.join([str(x) for x in row])+" " +
                mismatch_prim+" "+mismatch_hum+"\n")
    f.close()
