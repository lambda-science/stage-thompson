# Base de données Uniprot
cd /home/meyer/stage-thompson/data/raw/blastp-parameter-estimation/
database="callithrix_jacchus_uni.fasta  chlorocebus_sabaeus_uni.fasta macaca_fascicularis_uni.fasta nomascus_leucogenys_uni.fasta  pan_troglodytes_uni.fasta pongo_abelii_uni.fasta gorilla_gorilla_uni.fasta macaca_mulatta_uni.fasta papio_anubis_uni.fasta otolemur_garnettii_uni.fasta"
seqvar=$( cat *.id.results )
for db in $database
do
    echo db
    seqnum=$(echo $seqvar | grep -c "$db")
    echo "scale=1;$seqnum/500*100" | bc
done