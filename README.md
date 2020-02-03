# Code source pour mon stage de fin d'étude M2 Biotech HT Analyse @ ESBS Strasbourg - dans avec Julie D. Thompson dans l'équipe CSTB - iCube

### Structure des données (dossier local et sur serveur ssh.lbgi.fr)
```
├── bin  --------------------------------# Ensemble de script bash
├── doc  --------------------------------# Tout documents en rapport avec le projet (article, cahier de labo, notes)
├── raw  --------------------------------# Donnée brute ou en partie traitée
│   ├── orthoinspector-json -------------# Fichier raw JSON Orthoinspector
│   ├── orthoinspector-json-processed ---# Fichier issue du traitement des JSON
│   ├── uniprot-sequence ----------------# ID, Séquence, Alignement des protéines orthologue Uniprot
│   ├── uniprot-exon-map ----------------# Exon-map et autre séquence genomique, CDS, cDNA
│   ├── uniprot-error-mismatch ----------# Fichiers relatifs aux détection d'erreur d'alignement et mismatch
│   └── uniprot-translation-correction --# Fichiers relatifs à la correction des erreurs de mismatch par traduction
├── results  ----------------------------# Dossier contenant certains résultats finaux
├── temp  -------------------------------# Dossier pour tout fichier temporaire (screenshot, sequence test)
└── src  --------------------------------# Code source non-bash (python, jupyter notebook)
```
### Workflow pour re-générer les données
1. Executer src/1-Orthoinspector_to_ID_file.ipynb  
2. Executer bin/retrieve_seq_from_uniprot.sh  
3. Executer bin/mafft_all.sh   
4. Executer bin/error_count.sh
5. Executer ```cat ~/stage-thompson/raw/uniprot-exon-map/uniprot_errors.txt | grep "SEQ_ERROR3" > ~/stage-thompson/raw/uniprot-exon-map/uniprot_errors_mismatch.txt```
6. Récupérer que les ID des séquences à mismatch:  
```cat ~/stage-thompson/raw/uniprot-exon-map/uniprot_errors_mismatch.txt | cut -d ' ' -f3 | uniq > mismatch.id```
7. Récupérer les séquences des ID à mismatch
```/biolo/blast/bin/blastdbcmd -entry_batch mismatch.id -db /commun/bics/DB-Corentin/uniprot >> mismatch.id.fasta```
7. Optionnel: Executer src/2-Uniprot_RefSeq_Match.ipynb
8. Faire la correspondance UniprotID -> Ensembl Transcript avec https://www.uniprot.org/uploadlists/ (télécharger le fichier .tab)
9. Executer src/3-Retrieve_Genomic_CDS_and_Exon.ipynb
10. Executer src/4-Generate_Exon_Map.ipynb  
11. Executer src/5-ExonMap_to_cDNA.ipynb
12. Executer src/6-Mismatch_correction_translation.ipynb
