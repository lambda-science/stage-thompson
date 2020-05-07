FILES=~/stage-thompson/data/correction/two-by-two2/*.mafft
for f in $FILES
do
	/home/julie/ALN_UTILS/seqerrs $f >> ~/stage-thompson/data/correction/error_after_correction2.txt
done

FILES=~/stage-thompson/data/correction/two-by-two-old2/*.mafft
for f in $FILES
do
	/home/julie/ALN_UTILS/seqerrs $f >> ~/stage-thompson/data/correction/error_before_correction2.txt
done
