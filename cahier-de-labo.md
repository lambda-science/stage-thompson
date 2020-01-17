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
/biolo/blast/bin/blastp -db blast-gorilla -outfmt '6 sseqid sseq' -max_target_seqs 1 -query fakeseq.fasta | awk 'BEGIN{FS="\t"; OFS="\n"}{gsub(/-/, "", $2); print ">"$1,$2}'
```
