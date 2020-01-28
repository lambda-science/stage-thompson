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
