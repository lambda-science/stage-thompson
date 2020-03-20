FILES=~/stage-thompson/data/raw/uniprot-blast/*.mafft
for f in $FILES
do
	/home/julie/ALN_UTILS/seqerrs $f >> ~/stage-thompson/data/raw/uniprot_blast_errors.txt
done
