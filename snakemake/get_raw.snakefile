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

# Fichiers en sortie du workflow
# rule target:
#    input:
#        "../data/deletion-analysis/correction-pairwise/translation_match.tab",
#        "../data/deletion-analysis/uniprot-translation-correction/translation_match.tab",

rule target:
    input:
        "../data/raw/orthoinspector-json/*.json",
        "../data/raw/orthoinspector-json-processed/*.txt",
        "../data/raw/uniprot-sequence/*.id",
        "../data/raw/uniprot-sequence/*.id.fasta",
        "../data/raw/uniprot-sequence/*.id.fasta.mafft",
#############################################################################################

# 1/ Récupération des erreurs TYPE 2
rule get_errors:
    input:
        
    output:
        
    message:
        
    shell:
        
