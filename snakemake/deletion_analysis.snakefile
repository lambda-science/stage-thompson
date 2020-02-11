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
rule target:
    input:
        "../data/deletion-analysis/translation_match_100percent.tab",
        "../data/deletion-analysis/translation_match_similarity.tab",
        "../data/deletion-analysis/Exon_map.tab",
        "../data/deletion-analysis/Intron_map.tab"
#############################################################################################
# 9/ Correction erreur type2 par traduction et 80% similarité
rule error_correction_similarity_translation:
    input:
        "../data/deletion-analysis/uniprot_errors_type2.txt", "../data/deletion-analysis/genomic_all.fasta"
    params:
        out_folder="../data/deletion-analysis/similarity_correction/",
        mafft_path="/biolo/mafft/inst/bin/mafft"
        relatif="alignement"
    output:
        "../data/deletion-analysis/translation_match_similarity.tab"
    message:
        "Traduction des séquences et recherche du peptide human avec 80% similarité et <10 gap"
    shell:
        "python ../src/8-Mismatch_correction_similarity.py {input[0]} {input[1]} {params.out_folder} {params.mafft_path} {output} {params.relatif}"

rule error_correction_identity_translation:
    input:
        "../data/deletion-analysis/uniprot_errors_type2.txt", "../data/deletion-analysis/genomic_all.fasta"
    params:
        relatif="alignement"
    output:
        "../data/deletion-analysis/translation_match_100percent.tab"
    message:
        "Traduction des séquences génomique et recherche du peptide humain pour correction"
    shell:
        "python ../src/6-Mismatch_correction_translation.py {input[0]} {input[1]} {output} {params.relatif}"


#  5/ Création de l'exon_map
rule create_exon_map:
    input:
        "../data/deletion-analysis/transcript_ensembl.tab", "../data/deletion-analysis/exon_json_dump.json", "../data/deletion-analysis/genomic_all.fasta"
    output:
        "../data/deletion-analysis/Exon_map.tab", "../data/deletion-analysis/Intron_map.tab"
    message:
        "Generation de l'exon et de l'intron map"
    shell:
        "python ../src/4-Generate_Exon_Map.py {input[0]} {input[1]} {input[2]} {output[0]} {output[1]}"

#  4/ Récupération séquences génomique, CDS, info exon
rule get_genomic_CDS_exon_info:
    input:
        "../data/deletion-analysis/transcript_ensembl.tab"
    output:
        "../data/deletion-analysis/CDS_all.fasta", "../data/deletion-analysis/genomic_all.fasta", "../data/deletion-analysis/exon_json_dump.json"
    message:
        "Récupération seq génomique, seq CDS et json dump des exons"
    shell:
        "python ../src/3-Retrieve_genomic_CDS_and_Exon.pyled-1.py {input[0]} {output[0]} {output[1]} {output[2]}"

# 2.1/ Récupérer les séquences des ID à erreur
# 3/ Correspondance Uniprot ID -> Ensembl Transcript
rule get_SEQ_and_ID_errors:
    input:
        "../data/deletion-analysis/deletion.id"
    output:
        "../data/deletion-analysis/transcript_ensembl.tab", "../data/deletion-analysis/deletion.id.fasta"
    message:
        "Récupération des séquences correspondantes aux ID à erreurs\nRécupération identifiant Trasncript Ensembl pour chaque ID uniprot"
    shell:
        "python ../src/1.1-Uniprot-To-Ensembl.py {input} > {output[0]} &"
        "/biolo/blast/bin/blastdbcmd -entry_batch {input} -db /commun/bics/DB-Corentin/uniprot >> {output[1]}"

# 2/ Récupération des ID des erreurs
rule get_ID_errors:
    input:
        "../data/deletion-analysis/uniprot_errors_type2.txt"
    output:
        "../data/deletion-analysis/deletion.id"
    message:
        "Récupération des ID des erreurs"
    shell:
        "cat {input} | cut -d ' ' -f3 | uniq > {output}"

# 1/ Récupération des erreurs TYPE 2 et selection par taille (>=10)
rule get_errors:
    input:
        "../data/raw/uniprot_new_errors.txt"
    output:
        "../data/deletion-analysis/uniprot_errors_type2.txt"
    message:
        "Récupération des erreurs et selection par taille"
    shell:
        "cat {input} | grep \"SEQ_ERROR2\" | awk '$5-$4 >= 10 {print}' > {output}"
