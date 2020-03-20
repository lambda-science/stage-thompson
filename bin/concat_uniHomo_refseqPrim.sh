cd ~/stage-thompson/data/raw/refseq-blast/
FILES=*.fasta
for f in $FILES
do
	cat $f results/$f.results.fasta > results/$f.results.fasta2
done