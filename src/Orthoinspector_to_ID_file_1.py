# ## Script permettant de récupérer les Json Orthoinspector et d'obtenir un fichier par protéine humaine avec les identifiant Uniprot des protéines ortologue pour chaque protéine (humaine et primate)
# Input: aucun
# Output:
#     Stockage des réponse Json Orthoinspector
#     Stockage du processing des json avec un fichie par espèce contenant les protéines ID
#     Stockage d'un unique tableau contenant en colonne chaque espèce et en ligne chaque protéine ID
#     Stockage d'un fichier par ligne contenant l'ID humain et les ID de chaque primate

import pandas as pd
from os import listdir
from os.path import isfile, join
import json


def process_Json_orthoinspector():
    # Function: Process orthoinspector json response to make simple .tab files
    # Parameters: None
    # Return: None. Write files in raw/orthoinspector-json/
    # Description: Process orthoinspector json response to make simple .tab files with one file per primate vs human with a column
    # for human uniprot ID and one column for primate uniprot ID
    taxonID = [60711, 9541, 9544, 9555, 9595, 9598, 9601, 61853, 9483, 30611]
    for i in range(10):
        f = open("../data/raw/orthoinspector-json-processed-v2/" +
                 str(taxonID[i])+".txt", "w")
        with open("../data/raw/orthoinspector-json-v2/"+str(taxonID[i])+".json") as json_file:
            data = json.load(json_file)
            f.write("9606\t"+str(taxonID[i])+"\n")
            for human_prot in data:
                if human_prot == data[human_prot]['orthologs'][0]:
                    f.write(human_prot+"\t" + str("NaN"+"\n"))
                else:
                    f.write(human_prot+"\t" +
                            str(data[human_prot]['orthologs'][0]+"\n"))
        f.close()


def merge_tab_files():
    # Function: Merge processed .tab files from orthoinspector to form one single big table with all primates
    # Parameters: None
    # Return: main_df (dataframe): the pandas dataframe containing the table with all primates ID. Also write it to a file.
    # Description: Merged final final contains one column for each primate with each uniprot ID. Column 0 is always human.
    myfiles2 = listdir("../data/raw/orthoinspector-json-processed-v2/")
    main_df = pd.read_csv(
        '../data/raw/orthoinspector-json-processed-v2/'+myfiles2[0], sep='\t', header=0)
    for i in range(1, 10):
        df_temp = pd.read_csv(
            '../data/raw/orthoinspector-json-processed-v2/'+myfiles2[i], sep='\t', header=0)
        main_df = main_df.merge(df_temp, on="9606", how="outer")

    main_df.to_csv(
        "../data/raw/orthoinspector-json-processed-v2/allProteineName.txt", sep="\t", index=False)
    return main_df


def write_orthologue_ID(gene_names):
    # Function: Separate row of orthologous uniprot ID of all primate table to one row per file containing each ID.
    # Parameters:
    # 		gene_names: (dataframe) dataframe containing uniprotID for each human prot per row and each primates by columns
    # Return: None. Write a file for each row in /raw/uniprot-sequence/
    # Description:
    for i in range(20265):
        f = open("../data/raw/uniprot-sequence-v2/" +
                 str(gene_names.iloc[i, 0])+".id", "w", newline='\n')
        query = gene_names.iloc[i, :].dropna().tolist()
        query = ' '.join(query)
        query = query.split(" ")
        for i in query:
            f.write(str(i)+"\n")
        f.close()


if __name__ == "__main__":
    process_Json_orthoinspector()
    gene_names = merge_tab_files()
    write_orthologue_ID(gene_names)
