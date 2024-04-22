# make the G. hirsutum and G. barbadense transcriptome references from the G. raimondii reference genome (Paterson 2012)
curl -O https://www.cottongen.org/cottongen_downloads/Gossypium_raimondii/JGI_221_G.raimondii_Dgenome/assembly/G.raimondii_JGI_221_v2.0.assembly.fasta.gz
curl -O https://www.cottongen.org/cottongen_downloads/Gossypium_raimondii/JGI_221_G.raimondii_Dgenome/genes/G.raimondii_JGI_221_v2.1.transcripts.gff3.gz
gunzip *

sed -e '/>scaffold_14/,$d' -e 's/ ID=Chr.*$//g' G.raimondii_JGI_221_v2.0.assembly.fasta > Dgenome2_13.fasta
sed -e '/Gorai[.]N/d' G.raimondii_JGI_221_v2.1.transcripts.gff3 > D5.gff



# First, pseudogenomes of A2, AD1 At and Dt are constructed using `pseudogenome_by_snp.py`
module load py-biopython
python pseudogenome_by_snp.py Dgenome2_13.fasta snp41 pseudo4.1
python pseudogenome_by_snp.py Dgenome2_13.fasta snp42 pseudo4.2

# convert gff3 to gtf
module load cufflinks
gffread D5.gff -T -o D5.gtf


# append subgenome tag to Chrs and Gorai IDs
sed -e "s/-JGI_221_v2.1/.A/g" -e 's/\tphy/_A\tphy/' D5.gtf | grep "[.]A" > At.gtf
sed -e "s/-JGI_221_v2.1/.D/g" -e 's/\tphy/_D\tphy/' D5.gtf | grep "[.]D" > Dt.gtf

cat At.gtf Dt.gtf > pseudoAD.gtf

# built transcript reference from genome for RSEM
module load rsem
rsem-prepare-reference -p 20 --gtf pseudoAD.gtf pseudo4.1.At.fasta,pseudo4.1.Dt.fasta AD1_AtDt
rsem-prepare-reference -p 20 --gtf pseudoAD.gtf pseudo4.2.At.fasta,pseudo4.2.Dt.fasta AD2_AtDt

