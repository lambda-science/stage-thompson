# Code source pour mon stage de fin d'étude M2 Biotech HT Analyse @ ESBS Strasbourg - dans avec Julie D. Thompson dans l'équipe CSTB - iCube

### Workflow pour re-générer les données
1. Executer src/Orthoinspector_to_ID_file.ipynb  
2. Executer bin/retrieve_seq_from_uniprot.sh  
3. Executer bin/mafft_all.sh   
4. Executer bin/error_count.sh
5. Executer ```cat ~/stage-thompson/raw/uniprot-exon-map/uniprot_errors.txt | grep "SEQ_ERROR3" > ~/stage-thompson/raw/uniprot-exon-map/uniprot_errors_mismatch.txt```
6. Récupérer que les ID deséquences à mismatch:  
```cat ~/stage-thompson/raw/uniprot-exon-map/uniprot_errors_mismatch.txt | cut -d ' ' -f3 | uniq > mismatch.id```
7. Récupérer les séquences des ID à mismatch
```/biolo/blast/bin/blastdbcmd -entry_batch mismatch.id -db /commun/bics/DB-Corentin/uniprot >> mismatch.id.fasta```
