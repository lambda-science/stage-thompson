# Cahier de Laboratoire - suivis avancement du projet
***

### Semaine 1: 17.01.2020
* Récupération des identifiants uniprot des protéines orthologues homme / primates sur orthoinspector (script pyhon notebook orthologue.ipynb)
* Récupération des séquences de chaque protéines
* Alignement des 20000 fichiers de 11 séquences comparées avec MAFFT
* Bugfix variés (ex: Orthoinspectors qui indique la protéines recherché comme orthologue quand aucun orthologue n'est trouvé: source de répétitions de séquences)
* Familirisation avec le travail sur serveur distant ssh et utilisation de mini-script bash
* Projet: Amélioration des séquences en utilisant une base de données locale de uniprot récente
* Projet: Modifier la récupération des séquences en les récupérant grace à des blastP sur les protéomes des primates non-homme 
* Projet: Améliorer / Nettoyer l'organisation des fichiers sur ena/ et local.

Commande à retenir pour le blastp
```bash
/biolo/blast/bin/blastp -db blast-gorilla -evalue 0.005 -outfmt '6 sseqid sseq' -max_target_seqs 1 -query fakeseq.fasta | awk 'BEGIN{FS="\t"; OFS="\n"}{gsub(/-/, "", $2); print ">"$1,$2}'
```

* Projet: Récupérer BDD RefSeq. Faire des recherche gene name Humain - Primate OU/ET blast protéines humaine vers le best hit par primate.
* Récupérer exon map (Combien d'exon dans la région mal prédiet ? Pic du nb d'exon dans une window de 10 AA?) + carte genomique (Y'a des N dans les région mal prédite ?)
* A mettre à jour: récupérer blast-sw blast-trembl et re-faire les blastdbcmd pour récupérer les séquences à jour.
* A mettre à jour: refaire les requête orthoinspector pour récupérer les VRAI ortholog (bugfix signalé entrain d'être corrigé)

### Semaine 2: 20.01.2020
* Fait: installation des BDD RefSeq, Uniprot-SW, Uniprot-Trembl
* Fait: query des indentifiants uniprot sur la dernière version
* Fait: Nouveaux derniers alignements
* Projet: conversion Uniprot AC -> RefSeq Protein (beaucoup de perte)
* Projet: script conversion Uniprot AC -> Gene name -> RefSeq Protein (peu de perte mais très long et isoforme ?)
* Utilisation script julie compte d'erreur
* Projet: Faire la distribution du nb d'erreur par sequence -> normale ?
* Linéarisation de fasta et recherche d'un orga 
```cat * | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |grep "Macaca mu" -F -A 0 | tr "\t" "\n"```
* Annulation de cross-ref / blast chez RefSeq
* Simplement recherche des séquences présentant un mismatch chez Uniprot dans la base RefSeq avec 100% identité. Sur 30 000 seq avec mismatch uniprot, 5500 sont retrouvés avec 100% d'identité dans RefSeq. = Problême commun aux deux bases
* Erreur 3 = pas forcément une erreure de prédiction, peut être isoforme différent ? A confirmer par autre indices (NNN, exon ...). Regarder nb de iso ou de seq ds ref vs uni (= comparer qui en predit cb en moyenne). Missmatch = ok chez prot membranaire normalement (identité = trop strict). 30 000 seq = cb de fichie d'alignement ???
* 30 000 - 5500 de refgene = des isoformes unique à Uniprot
* Projet: récupérer les séquences de transcript (uniprot mapping) et récupérer donc les Exons. Pour les zones à erreur n°3 regarder combien d'exon son colocalisé (+2 ou 0 = erreur de prédiction) et regarder si présence de N dans seq génomique colocalisée avec erreur.

### Semaine 3: 27.01.2020
* En cours: requête des séquences génomiques des protéines présentant une erreur de type 3 (30 000 séquences à récupérer, 600 POST de 50 seq avec 30 en une fois)
* Si Error type 3: voir si on peut retrouver la séquence protéique humaine dans la zone génomique primate ça veut dire pb de prediction (Difficile si error 3 en N-ter mais faisable si entre deux exons)
* Modif sensibilité mismatch taille mini 10 AA -> 11659 seq
* Done: récupération exon map, seq genomique, cDNA. En cours: CDS
* Done: compter les petits exon et les seq cDNA contenant des N
* A faire: intron map, grep motif intron canonique cb ne le respecte pas ?
* A faire: retrouver seq protéique humaine dans seq génomique ? -> Algorithmie pour lire mafft, fichier erreur ect...
* Compter les séquences comportant de N 
```
cat cds_new.fa | grep -E "^[^>].*" | grep -E "N" | wc -l
```

**Discussion julie 29/01:**

* TODO: colocaliser les petits exons avec l'erreur
* TODO: compter les séquences ayant des XXXXX au mismatch
* TODO: traduction de zone prot humaine -> génomique
* TODO: Comparaison  stat exon/N avec gènes sans erreur/sans mismatch / randomly sampled
* J'ai vu des plusieurs petits exon séparés par 1 ou 2 bp introns -> intéressant détections d'introns sous la taille connue ?
* Arabidopsis taliana: introns 59 bp min, exon 1 bp min
* Humain mini 39 bp introns ? 
* Big TODO: refactorer et commenter le code pour le ré-utiliser
* Refaire la pipeline en filtrant les données pour les seq humaine n'existant pas ?
* Prototype de traduction: faire alignement MAAFT + score ?
* Journée du 30/01: grand netoyage du dossier local et SSH: fusion des deux et unificiation, supression des dossier et fichier inutile (création d'un back-up avant). Push et création d'un dépot git avec les sources. Organisation des scripts python, ajout de quelques commentaire et du détail du workflow necessaire pour re-creer les données. Ajout de doc sur la structure du dossier-projet.

### Semaine 4: 03/02/2020
* Update git pour utiliser des EOL LF seulement et ne plus commit les file permissions THANKS GOD  
* Réu julie/nico: discussion sur seq retrouvée. -> bug du a un frameshift de 1 NT, exon bon mais prot décalée ??? etrange + 1 exon en plus pourquoi ? Pas de grande dif de seq ou de splicing site  
* Fait: script traduction, correction exon par similarité  
* TO-DO: Regarder quel primate présente le plus d'erreur  
* Refaire l'analyse pour les erreurs de type 2 (deletion en milieu de séquence = exon loupé ?) -> traduction et localisation exon manquant  
* Regarder +/- 2 nt autour des exons prédit de correction pour voir si présence d'introns non-canonique.
* Done: error par primate
* En cours: projet snakemake 
* En cours: documentation du code  
* En cours: analyse erreur type 2

### Semaine 5: 10/02/2020
* Done: snakemake pour erreur type 2 (deletion)
* Done: doc de la moitié du code
* Done: Diapo présentation groupe A-M
* En cours: snakemake qui tourne
* En cours: préparation snakemake pour données raw
* Olivier: regarder seq à erreur uniprot spécific ET ou l'on retrouve pas la bonne séquence (gap) -> est-ce que le NCBI prédit aussi des trucs bidons ? -> DONE
* Done: Biblio sur Genemark Hmm ES 3.0  
* A faire: wildcare snakemake raw
* BUG FIX: Erreur Type 2: les pos sont des positions sur l'alignement et non par sur la seq humaine
* Réunion avec Julie parler du plan d'explications des mismatch. Info sur repeats, isoforme, seq avec N.  
* A programmer: Flagging des mismatch, détection des 7 features différentes (liste carnet)  

### Semaine 6: 17/02/2020
* Done: 4/7 flagging function
* In progress: selecting the correct ensembl transcript ID 
* Journal club article
* In progress: mismatch exon map for human  
* Chemin fichier kirsley /gstock/user/kchennen/Pampas/data/1_interim/Homo_sapiens/primates/
