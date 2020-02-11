########################################################
####### Récupération des données brute l'analyse #######
########################################################
# input: None.
#
# output:
#        "../data/raw/orthoinspector-json/*.json",
#        "../data/raw/orthoinspector-json-processed/*.txt",
#        "../data/raw/uniprot-sequence/*.id",
#        "../data/raw/uniprot-sequence/*.id.fasta",
#        "../data/raw/uniprot-sequence/*.id.fasta.mafft",
#        "../data/raw/uniprot-sequence/uniprot_new_errors.txt"
#
# processus :
#       1/ Récupération données orthoinspectors
#       2/ Processing données orthoinspectors
#       3/ Récupération des séquences.
#       4/ Alignement des séquences.
#       5/ Error-calling sur les alignements.

#################################################
#################################################

rule target:
    input:
        "../data/raw/orthoinspector-json/",
        "../data/raw/orthoinspector-json-processed/",
        "../data/raw/uniprot-sequence/*.id",
        "../data/raw/uniprot-sequence/*.id.fasta",
        "../data/raw/uniprot-sequence/*.id.fasta.mafft",
        "../data/raw/uniprot-sequence/uniprot_new_errors.txt"

#############################################################################################
# 5/ Error-calling sur les alignements.
rule error_calling_julies_script: 
    input:
        "../data/raw/uniprot-sequence/*.id.fasta.mafft"       
    output: 
        "../data/raw/uniprot-sequence/uniprot_new_errors.txt"   
    message: 
        "Détermination des erreurs d'alignement"
    shell: 
        "./../bin/orthoinspect_query.sh"
        
# 4/ Alignement des séquences.
rule mafft_align:
    input:
        "../data/raw/uniprot-sequence/*.id.fasta"        
    output: 
        "../data/raw/uniprot-sequence/*.id.fasta.mafft"
    message: 
        "Alignement des séquences orthologues"
    shell: 
        "./../bin/mafft_all.sh"
        
# 3/ Récupération des séquences.
rule DB_get_seq:
    input: 
        "../data/raw/uniprot-sequence/*.id"    
    output: 
        "../data/raw/uniprot-sequence/*.id.fasta"
    message: 
        "Récupération séquences protéïque dans la BDD Uniprot"
    shell: 
        "./../bin/retrieve_seq_from_uniprot.sh"
        
# 2/ Processing données orthoinspectors
rule orthoinspector_process:
    input:
        "../data/raw/orthoinspector-json/"        
    output: 
        directory("../data/raw/orthoinspector-json-processed"),
        "../data/raw/uniprot-sequence/*.id"
    message: 
        "Processing des données Orthoinspector"
    shell: 
        "python ../src/1-Orthoinspector_to_ID_file.py"
        
# 1/ Récupération données orthoinspectors
rule orthoinspector_requests:        
    output: 
        directory("../data/raw/orthoinspector-json") 
    message: 
        "Récupération des données Orthoinspector"
    shell: 
        "./../bin/orthoinspect_query.sh"
        