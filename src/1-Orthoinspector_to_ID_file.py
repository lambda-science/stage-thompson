# ## Script permettant de récupérer les Json Orthoinspector et d'obtenir un fichier par protéine humaine avec les identifiant Uniprot des protéines ortologue pour chaque protéine (humaine et primate)
# Input: aucun  
# Output:  
#     Stockage des réponse Json Orthoinspector  
#     Stockage du processing des json avec un fichie par espèce contenant les protéines ID  
#     Stockage d'un unique tableau contenant en colonne chaque espèce et en ligne chaque protéine ID  
#     Stockage d'un fichier par ligne contenant l'ID humain et les ID de chaque primate

# Traiter les fichiers json pour avoir un colonne prot humaine et une colonne orthologue chez espèce d'intérêt
import pandas as pd
from os import listdir
from os.path import isfile, join
import json

def process_Json_orthoinspector():
    # Traiter les réponses json orthoinspector pour en faire des fichier tabulés simple
    taxonID = [60711, 9541, 9544, 9555, 9595, 9598, 9601, 61853, 9483, 30611]
    human = 9606

    myfiles = listdir("../raw/orthoinspector-json/")
    for i in range(10):
        f = open("../raw/orthoinspector-json-processed/"+str(taxonID[i])+".txt", "w")
        with open("../raw/orthoinspector-json/"+str(taxonID[i])+".json") as json_file:
            data = json.load(json_file)
            f.write("9606\t"+str(taxonID[i])+"\n")
            for human_prot in data:
                if human_prot == data[human_prot]['orthologs'][0]:
                    f.write(human_prot+"\t"+ str("NaN"+"\n"))
                else:
                    f.write(human_prot+"\t"+ str(data[human_prot]['orthologs'][0]+"\n"))
        f.close()

def merge_tab_files():
    # Fusion des fichier traiter précédent pour passer d'un fichier par comparaison à un fichier pour toute les comparaison
    myfiles2 = listdir("../raw/orthoinspector-json-processed/")
    main_df = pd.read_csv('../raw/orthoinspector-json-processed/'+myfiles2[0], sep='\t', header=0)
    for i in range(1, 10):
        df_temp = pd.read_csv('../raw/orthoinspector-json-processed/'+myfiles2[i], sep='\t', header=0)
        main_df = main_df.merge(df_temp, on="9606", how="outer")
    
    main_df.to_csv("../raw/orthoinspector-json-processed/allProteineName.txt", sep="\t", index=False)
    return main_df

def write_orthologue_ID(gene_names):
    # Préparation des identifiants pour query de séquence par fichier. Un fichier par comparaison (20 000 au total) à récupérer
    # Une séquence par ligne
    for i in range(20265):
        f = open("../raw/uniprot-sequence/"+str(gene_names.iloc[i,0])+".id", "w", newline='\n')
        query = gene_names.iloc[i,:].dropna().tolist()
        query = ' '.join(query)
        query = query.split(" ")
        for i in query:
            f.write(str(i)+"\n")
        f.close()

if __name__=="__main__":
    process_Json_orthoinspector()
    gene_names = merge_tab_files()
    write_orthologue_ID(gene_names)
