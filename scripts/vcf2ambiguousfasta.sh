#!/bin/bash
set -e # makes errors loud and noticeable
set -u # aborting the script if a variable’s value is unset
set -o pipefail # makes errors in a pipe loud and noticeable

# Script to cut genes out of a vcf and create fasta files for each them using ambiguity codes.

# create folders
# "originals" contains all the raw input files
# "tmp" contains modified input files
# "fasta_ingroup" contains the fasta files for each gene of all samples of Oc (without reference)
# "fasta_reference" contains the fasta files for each gene of reference and outgroup
mkdir originals
mkdir tmp
mkdir fasta_ingroup
mkdir -p fasta_references/{fasta_FIOER33,fasta_BEOM2}

# copy all needed data to originals folder:
# reference and annotation file (FIOER33 & BEOM2)
cp /home/pascal/master_publication/data/FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta originals/.
cp /home/pascal/master_publication/data/FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.gff3 originals/.
cp /home/pascal/master_thesis/data/reference_genome_NCBI/BEOM2/GCA_004325065.1_BEOM2_v1_cds_from_genomic.fna originals/.
# perl script to cut the gene sequences out of the genome-fasta according to the gff
cp /home/pascal/master_publication/scripts/gff2fasta_Ham.pl originals/.
# vcf

# originals from orthomcl

#!! define samples to use !!
samples=$(bcftools query -l ../hamiltosporidium_10122020_SNP.recode.vcf) # all samples in the vcf

# prepare the reference genome for gatk
picard CreateSequenceDictionary R=originals/FI-OER-3-3.pacbio.v1.fasta O=originals/FI-OER-3-3.pacbio.v1.dict
samtools faidx originals/FI-OER-3-3.pacbio.v1.fasta

##########################################################################################
# Samples #
###########

# loop to generate the fasta files of all samples in the vcf for each gene
for sample in $samples
    do

    # new vcf containing only one sample
    GenomeAnalysisTK -T SelectVariants -R originals/FI-OER-3-3.pacbio.v1.fasta -o $sample.vcf -V ../hamiltosporidium_10122020_SNP.recode.vcf -sn $sample
    
    #
    sed -i "" '/*/d' $sample.vcf

    # new "reference" for the one sample (with ambiguities)
    GenomeAnalysisTK -T FastaAlternateReferenceMaker -R originals/FI-OER-3-3.pacbio.v1.fasta -o $sample.fasta -V $sample.vcf --use_IUPAC_sample $sample

    # Replace blanks by dashes. Works for: up to three gaps in a row
    sed -E -e "s/   ([^c ])/---\1/g;s/  ([^c ])/--\1/g;s/ ([^c ])/-\1/g;s/([A-Z])   $/\1---/g;s/([A-Z])  $/\1--/g;s/([A-Z]) $/\1-/g" $sample.fasta > $sample.gap-.fasta

	sed -i "" -E "s/>.*-(.*):1/>\1/g" $sample.gap-.fasta

    # cut the gene sequences out of the genome-fasta according to the gff for the sample
    perl originals/gff2fasta_Ham.pl $sample.gap-.fasta originals/FI-OER-3-3.pacbio.v1.gff3 $sample

    # prepare new sample fasta (only containing coding sequences)
    sed -i "" 's/gnl|WGS:PITJ|//g' $sample.cds.fasta
    sed -i "" -E "s/(>.*)/\1|$sample/g" $sample.cds.fasta
    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $sample.cds.fasta  > $sample.cds.oneline.fasta
    
    # cut orthologs
    for i in `cat ../../../cds/orthomcl/orthomcl-output_25012021/groups/orthologs_of_FIOER33_to_BEOM2.txt`;do awk 'BEGIN{RS=">"} /'$i'/{print ">" $0}' $sample.cds.oneline.fasta;done > $sample.orthologs.fasta

    # create a fasta file for each of the genes
    awk -F "|" '/^>/ {close(F) ; F = $1".fasta"} {print >> F}' $sample.orthologs.fasta

    # move all unneeded files to the tmp folder
    mv $sample* tmp/

done

# move the gene-fastas to the fasta_ingroup folder
for f in *FUN*; do mv "$f" "$(echo "fasta_ingroup/$f" | sed -e 's/>//')";done

##########################################################################################
# Reference plus outgroup #
###########################

# cut the gene sequences out of the reference according to the gff
perl originals/gff2fasta_Ham.pl originals/FI-OER-3-3.pacbio.v1.fasta originals/FI-OER-3-3.pacbio.mod.gff3 Reference

# prepare reference fasta (only containing coding sequences)
sed -i "" -E "s/(>.*)/\1|FIOER33/g" Reference.cds.fasta
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Reference.cds.fasta  > Reference.cds.oneline.fasta

# cut orthologs
for i in `cat ../../../cds/orthomcl/orthomcl-output_25012021/groups/orthologs_of_FIOER33_to_BEOM2.txt`;do awk 'BEGIN{RS=">"} /'$i'/{print ">" $0}' Reference.cds.oneline.fasta;done > Reference.orthologs.fasta

# create a fasta file for each of the genes
awk -F "|" '/^>/ {close(F) ; F = $1".fasta"} {print >> F}' Reference.orthologs.fasta

# move the gene-fastas to the fasta_references folder
for f in *FUN*; do cp "$f" "$(echo "fasta_references/fasta_FIOER33/$f" | sed -e 's/>//')";done

# prepare the outgroup fasta (only containing coding sequences)
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < originals/GCA_004325065.1_BEOM2_v1_cds_from_genomic.fna > cds_BEOM2.oneline.fna
sed -i "" -E 's/(>)lcl.PI............_cds_(..........).*/\1\2|BEOM2/g' cds_BEOM2.oneline.fna

# cut orthologs
for i in `cat ../../../cds/orthomcl/orthomcl-output_25012021/groups/orthologs_of_BEOM2_to_FIOER33.txt`;do awk 'BEGIN{RS=">"} /'$i'/{print ">" $0}' cds_BEOM2.oneline.fna;done > BEOM2.orthologs.fasta

# translate gene IDs from BEOM2 to FIOER33
paste ../../../cds/orthomcl/orthomcl-output_25012021/groups/orthologs_of_BEOM2_to_FIOER33.txt ../../../cds/orthomcl/orthomcl-output_25012021/groups/orthologs_of_FIOER33_to_BEOM2.txt | while read a b; do sed -i "" "s/$a/$b/" BEOM2.orthologs.fasta; done

# add the genes to the fastas of the reference genes
awk -F "|" '/^>/ {close(F) ; F = $1".fasta"} {print >> F}' BEOM2.orthologs.fasta
for f in *FUN*; do mv "$f" "$(echo "fasta_references/$f" | sed -e 's/>//;s/.fasta//')";done # as input list for Ka/Ks calculation

# move the gene-fastas to the fasta_references folder and the rest to the tmp folder
awk -F "|" '/^>/ {close(F) ; F = $1".fasta"} {print >> F}' BEOM2.orthologs.fasta
for f in *FUN*; do mv "$f" "$(echo "fasta_references/fasta_BEOM2/$f" | sed -e 's/>//')";done

# move all unneeded files to the tmp folder
mv Reference* tmp/
mv *BEOM2* tmp/
mv *shortHeader* tmp/
