#################################################
####### Analyse des alignement à deletion #######
#################################################
# input:
#      dossier raw/
#
# output:
#       ../data/deletion-analysis/correction-pairwise/translation_match.tab
#       ../data/deletion-analysis/uniprot-translation-correction/translation_match.tab
#
# processus :
#     1/ Récupération des erreurs
#     2/ Récupération des ID des erreurs
#     2.1 / Récupérer les séquences des ID à mismatch
#     3/ Correspondance Uniprot ID -> Ensembl Transcript
#     4/ Récupération séquences génomique, CDS, info exon
#     5/ Création de l'exon_map
#     ((6/ Création cDNA))
#     7/ Correction mismatch 100% identitié
#     8/ Colocalization mismatch exon/intron
#     9/ Correction mismatch par similarité
#
#################################################
#################################################

rule target:
    input:
        "../data/raw/orthoinspector-json/*.json",
        "../data/raw/orthoinspector-json-processed/*.txt",
        "../data/raw/uniprot-sequence/*.id",
        "../data/raw/uniprot-sequence/*.id.fasta",
        "../data/raw/uniprot-sequence/*.id.fasta.mafft",
        "../data/raw/uniprot-sequence/uniprot_new_errors.txt"

#############################################################################################
# 5/ Error-calling sur les alignements.
rule error_calling_julies_script:        
    output: "../data/raw/uniprot-sequence/uniprot_new_errors.txt"   
    message: "Détermination des erreurs d'alignement"
    shell: "./../bin/orthoinspect_query.sh"
        
# 4/ Alignement des séquences.
rule mafft_align:        
    output: "../data/raw/uniprot-sequence/*.id.fasta.mafft"   
    message: "Alignement des séquences orthologues"
    shell: "./../bin/mafft_all.sh"
        
# 3/ Récupération des séquences.
rule DB_get_seq:        
    output: "../data/raw/uniprot-sequence/*.id.fasta"   
    message: "Récupération séquences protéïque dans la BDD Uniprot"
    shell: "./../bin/retrieve_seq_from_uniprot.sh"
        
# 2/ Processing données orthoinspectors
rule orthoinspector_process:        
    output: "../data/raw/orthoinspector-json-processed/*.txt", "../data/raw/uniprot-sequence/*.id.fasta"   
    message: "Processing des données Orthoinspector"
    shell: "python ../src/1-Orthoinspector_to_ID_file.py"
        
# 1/ Récupération données orthoinspectors
rule orthoinspector_requests:        
    output: "../data/raw/orthoinspector-json/*.json"   
    message: "Récupération des données Orthoinspector"
    shell: "./../bin/orthoinspect_query.sh"
        