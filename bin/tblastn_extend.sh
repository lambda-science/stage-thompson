while read p; do
  /biolo/blastplus/blast-2.9.0/bin/tblastn -query $1/tblastn2/query_subject/$p.query -subject $1/tblastn2/query_subject/$p.subject -evalue 0.005 -max_target_seqs 1 -outfmt '6 sseqid qseqid sseq qseq sstart send evalue positive length sframe' > $1/tblastn2/blast_out/$p.blast 2>/dev/null
done < $1/tblastn2/all_couple.txt
touch $1/tblastn2/tblastn.done
