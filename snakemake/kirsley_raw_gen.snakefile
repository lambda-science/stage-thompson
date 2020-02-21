#########################################################
####### Analyse des alignement à mismatch kirsley #######
#########################################################

rule target:
    input:
        "../data/raw/kirsley_error.txt",
        "../data/raw/kirsley_raw.done"

# 0.5/ Error-calling sur les alignements.
rule error_calling_julies_script: 
    input:
        "../data/raw/kirsley/"    
    output: 
        "../data/raw/kirsley_error.txt"   
    message: 
        "Détermination des erreurs d'alignement"
    shell: 
        "./../bin/error_count_kirsley.sh"

        
# 0/ Alignement des séquences.
rule mafft_align:
    input:
        "../data/raw/kirsley/"    
    output: 
        "../data/raw/kirsley/*.fasta.mafft"
    message: 
        "Alignement des séquences orthologues"
    shell: 
        "./../bin/mafft_kirsley.sh & touch ../data/raw/kirsley_raw.done"
