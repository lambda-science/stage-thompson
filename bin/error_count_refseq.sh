FILES=~/stage-thompson/data/raw/refseq-blast/*.mafft
for f in $FILES
do
	/home/julie/ALN_UTILS/seqerrs $f >> ~/stage-thompson/data/raw/refseq_error.txt
done
