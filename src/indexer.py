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
    dictionary = dict(zip(seq, title))
    return dictionary

def importQuery(pathFasta):
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

my_DB = fasta2List("/commun/bics/DB-Corentin/refseq-perfect/all.fasta")
my_Query = importQuery("/home/meyer/refseq-perfect-match/ID_error3_filt.fasta")
f = open("results.out","w")
for i, j in my_Query.items():
    try: 
        match = my_DB[j]
        f.write("Match found: "+i+" "+my_DB[j]+"\n")
    except:
        pass
f.close()
