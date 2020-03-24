FILES=/home/meyer/stage-thompson/data/raw/refseq-blast/results/*.fasta
for f in $FILES
do
    sed -i 's/>.*|\(.*\)|/>\1 /' $f
done
