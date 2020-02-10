FILES=~/stage-thompson/data/raw/uniprot-sequence/*.mafft
for f in $FILES
do
	/home/julie/ALN_UTILS/seqerrs $f >> ~/stage-thompson/data/raw/uniprot_errors.txt
done
