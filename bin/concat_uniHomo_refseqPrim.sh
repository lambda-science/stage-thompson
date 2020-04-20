cd ~/stage-thompson/data/raw/refseq-blast/
FILES=*.id
for f in $FILES
do
	cat $f results/$f.fasta > results/$f.fasta2
done