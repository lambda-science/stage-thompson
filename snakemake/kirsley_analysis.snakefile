#########################################################
####### Analyse des alignement à mismatch kirsley #######
#########################################################

rule target:
    input:
        "../data/raw/kirsley/*.id.fasta",
        "../data/raw/kirsley/*.id.fasta.mafft",
        "../data/raw/kirsley/kirsley_error.txt",
        "../data/kirsley-analysis/kirsley_error_type3.txt",
        "../data/kirsley-analysis/blast_out/",
        "../data/kirsley-analysis/query_subject/",
        "../data/kirsley-analysis/match.out,
        "../data/kirsley-analysis/all_couple.txt",
        "../data/kirsley-analysis/Exon_map.tab",
        "../data/kirsley-analysis/Intron_map.tab"

##############################################################################
# 7/ Colocalize
rule mismatch_localize_exon:
    input:
        ""        
    output: 
        ""
    message: 
        "Colocalisation des mismatch sur l'exon map"
    shell: 
        ""
# 6.2/ tBlastn process
rule tblastn_process:
    input:
        "../data/mismatch-analysis/tblastn/blast_out/"        
    output: 
        "../data/mismatch-analysis/tblastn/blast_out/"
    message: 
        "Process des résultats tblastn"
    shell: 
        "./../bin/tblastn_kirsley.sh"

# 6.1/ tBlastn exec
rule tblastn_exec:
    input:
        "../data/mismatch-analysis/tblastn/all_couple.txt"        
    output: 
        "../data/mismatch-analysis/tblastn/blast_out/"
    message: 
        "tblastn en cours"
    shell: 
        "./../bin/tblastn_kirsley.sh"

# 6/ tBlastn prepare
rule tblastn_prepare:
    input:
        "../data/mismatch-analysis/tblastn/all_couple.txt"        
    output: 
        "../data/mismatch-analysis/tblastn/blast_out/"
    message: 
        "préparation tblastn"
    shell: 
        "./../bin/tblastn_kirsley.sh"

# 5/ Création de l'exon_map
rule create_exon_map:
    input:
        "../data/kirsley-analysis/transcript_ensembl_corrected.tab", "../data/kirsley-analysis/exon_json_dump.json", "../data/kirsley-analysis/genomic_all.fasta"
    output:
        "../data/kirsley-analysis/Exon_map.tab", "../data/kirsley-analysis/Intron_map.tab"
    message:
        "Generation de l'exon et de l'intron map"
    shell:
        "python ../src/Generate_Exon_Map_4.py {input[0]} {input[1]} {input[2]} {output[0]} {output[1]}"

#  4/ Récupération séquences génomique, CDS, info exon
rule get_genomic_CDS_exon_info:
    input:
        "../data/kirsley-analysis/transcript_ensembl_corrected.tab"
    output:
        "../data/kirsley-analysis/CDS_all.fasta", "../data/kirsley-analysis/genomic_all.fasta", "../data/kirsley-analysis/exon_json_dump.json"
    message:
        "Récupération seq génomique, seq CDS et json dump des exons"
    shell:
        "python ../src/Retrieve_genomic_CDS_and_Exon_3.py {input[0]} {output[0]} {output[1]} {output[2]}"

# 2.1/ Récupérer les séquences des ID à erreur
# 3/ Correspondance Uniprot ID -> Ensembl Transcript
# ICI: grep dans fasta
# AJOUTER ETAPE SELECTION
rule get_SEQ_and_ID_errors:
    input:
        "../data/kirsley-analysis/mismatch.id"
    output:
        "../data/kirsley-analysis/transcript_ensembl.tab", "../data/kirsley-analysis/mismatch.id.fasta"
    message:
        "Récupération des séquences correspondantes aux ID à erreurs\nRécupération identifiant Trasncript Ensembl pour chaque ID uniprot"
    shell:
        "python ../src/Uniprot-To-Ensembl_1_1.py {input} > {output[0]} &"
        "/biolo/blast/bin/blastdbcmd -entry_batch {input} -db /commun/bics/DB-Corentin/uniprot >> {output[1]}"

# 2/ Récupération des ID des erreurs
rule get_ID_errors:
    input:
        "../data/kirsley-analysis/kirsley_error_type3.txt"
    output:
        "../data/kirsley-analysis/mismatch.id"
    message:
        "Récupération des ID des erreurs"
    shell:
        "cat {input} | cut -d ' ' -f3 | uniq > {output}"

# 1/ Récupération des erreurs TYPE 3
rule get_errors:
    input:
        "../data/raw/kirsley_error.txt"
    output:
        "../data/kirsley-analysis/kirsley_error_type3.txt"
    message:
        "Récupération des erreurs"
    shell:
        "cat {input} | grep \"SEQ_ERROR3\" > {output}"

# 0.5/ Error-calling sur les alignements.
rule error_calling_julies_script: 
    input:
        "../data/raw/kirsley/*.id.fasta.mafft"       
    output: 
        "../data/raw/kirsley/kirsley_error.txt"   
    message: 
        "Détermination des erreurs d'alignement"
    shell: 
        "./../bin/error_count_kirsley.sh"
        
# 0/ Alignement des séquences.
rule mafft_align:
    input:
        "../data/raw/kirsley/*.id.fasta"        
    output: 
        "../data/raw/kirsley/*.id.fasta.mafft"
    message: 
        "Alignement des séquences orthologues"
    shell: 
        "./../bin/mafft_kirsley.sh"
