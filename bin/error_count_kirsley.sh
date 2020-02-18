FILES=~/stage-thompson/data/raw/kirsley/*.mafft
for f in $FILES
do
	/home/julie/ALN_UTILS/seqerrs $f >> ~/stage-thompson/data/raw/kirsley_error.txt
done
