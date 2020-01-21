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

