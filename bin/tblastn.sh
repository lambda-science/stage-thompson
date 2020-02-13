while read p; do
  /biolo/blastplus/blast-2.9.0/bin/tblastn -query data/mismatch-analysis/tblastn/query_subject/$p.query -subject data/mismatch-analysis/tblastn/query_subject/$p.subject -evalue 0.005 -max_target_seqs 1 -outfmt '6 sseqid qseqid sseq qseq sstart send evalue positive length sframe' > data/mismatch-analysis/tblastn/blast_out/$p.blast
done < data/mismatch-analysis/tblastn/all_couple.txt
