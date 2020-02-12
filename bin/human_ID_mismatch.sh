grep "SEQ_ERROR3" data/raw/uniprot_new_errors.txt | grep ".*.id" -o | \
sed "s/\/home\/meyer\/uniprot\/\(.*\)\.id/\1/" | uniq \
> data/mismatch-flagging/human_uniprot_ID.id
python src/Uniprot_To_Ensembl_1_1.py data/mismatch-flagging/human_uniprot_ID.id > data/mismatch-flagging/human_uniprot_ensembl.tab
