# variant calling of 2 samples at a time
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta FI-BR1-39-6.0.23.mapped.md.bam FI-OER-3-3.0.05.mapped.md.bam | bcftools call -c - > pseudo-diploid/FI-BR1-39-6_FI-OER-3-3.bcftools.vcf
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta FI-BR1-39-6.0.23.mapped.md.bam RU-TY1-3.0.13.mapped.md.bam | bcftools call -c - > pseudo-diploid/FI-BR1-39-6_RU-TY1-3.bcftools.vcf
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta FI-BR1-39-6.0.23.mapped.md.bam IL-SH-4.0.08.mapped.md.bam | bcftools call -c - > pseudo-diploid/FI-BR1-39-6_IL-SH-4.bcftools.vcf
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta FI-BR1-39-6.0.23.mapped.md.bam ES-OM-5.0.14.mapped.md.bam | bcftools call -c - > pseudo-diploid/FI-BR1-39-6_ES-OM-5.bcftools.vcf
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta FI-BR1-39-6.0.23.mapped.md.bam BE-OM-2.0.08.mapped.md.bam | bcftools call -c - > pseudo-diploid/FI-BR1-39-6_BE-OM-2.bcftools.vcf
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta IL-SH-4.0.08.mapped.md.bam BE-OM-2.0.08.mapped.md.bam | bcftools call -c - > pseudo-diploid/IL-SH-4_BE-OM-2.bcftools.vcf
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta ES-DO1-1.0.18.mapped.md.bam BE-OM-2.0.08.mapped.md.bam | bcftools call -c - > pseudo-diploid/ES-DO1-1_BE-OM-2.bcftools.vcf
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta RU-BAYA1-1.0.07.mapped.md.bam BE-OM-2.0.08.mapped.md.bam | bcftools call -c - > pseudo-diploid/RU-BAYA1-1_BE-OM-2.bcftools.vcf
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta IL-SH-4.0.08.mapped.md.bam ES-DO1-1.0.18.mapped.md.bam | bcftools call -c - > pseudo-diploid/IL-SH-4_ES-DO1-1.bcftools.vcf
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta FI-BR1-39-6.0.23.mapped.md.bam RU-BAYA1-1.0.07.mapped.md.bam | bcftools call -c - > pseudo-diploid/FI-BR1-39-6_RU-BAYA1-1.bcftools.vcf
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta IL-SH-4.0.08.mapped.md.bam RU-BAYA1-1.0.07.mapped.md.bam | bcftools call -c - > pseudo-diploid/IL-SH-4_RU-BAYA1-1.bcftools.vcf
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta ES-DO1-1.0.18.mapped.md.bam RU-BAYA1-1.0.07.mapped.md.bam | bcftools call -c - > pseudo-diploid/ES-DO1-1_RU-BAYA1-1.bcftools.vcf
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta BE-OM-2.0.08.mapped.md.bam BE-T1-2.0.35.mapped.md.bam | bcftools call -c - > pseudo-diploid/BE-OM-2_BE-T1-2.bcftools.vcf
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta BE-OM-2.0.08.mapped.md.bam BE-WH2-11.0.34.mapped.md.bam | bcftools call -c - > pseudo-diploid/BE-OM-2_BE-WH2-11.bcftools.vcf
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta BE-T1-2.0.35.mapped.md.bam BE-WH2-11.0.34.mapped.md.bam | bcftools call -c - > pseudo-diploid/BE-WH2-11_BE-T1-2.bcftools.vcf

# index BAM files
for md in *md.bam; do samtools index $md; done

cd pseudo-diploid/

# read-based phasing
whatshap phase --reference ../../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta -o FI-BR1-39-6_FI-OER-3-3.phased.bcftools.vcf FI-BR1-39-6_FI-OER-3-3.bcftools.vcf ../FI-BR1-39-6.0.23.mapped.md.bam ../FI-OER-3-3.0.05.mapped.md.bam
whatshap phase --reference ../../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta -o FI-BR1-39-6_RU-TY1-3.phased.bcftools.vcf FI-BR1-39-6_RU-TY1-3.bcftools.vcf ../FI-BR1-39-6.0.23.mapped.md.bam ../RU-TY1-3.0.13.mapped.md.bam
whatshap phase --reference ../../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta -o FI-BR1-39-6_IL-SH-4.phased.bcftools.vcf FI-BR1-39-6_IL-SH-4.bcftools.vcf ../FI-BR1-39-6.0.23.mapped.md.bam ../IL-SH-4.0.08.mapped.md.bam
whatshap phase --reference ../../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta -o FI-BR1-39-6_ES-OM-5.phased.bcftools.vcf FI-BR1-39-6_ES-OM-5.bcftools.vcf ../FI-BR1-39-6.0.23.mapped.md.bam ../ES-OM-5.0.14.mapped.md.bam
whatshap phase --reference ../../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta -o FI-BR1-39-6_BE-OM-2.phased.bcftools.vcf FI-BR1-39-6_BE-OM-2.bcftools.vcf ../FI-BR1-39-6.0.23.mapped.md.bam ../BE-OM-2.0.08.mapped.md.bam
whatshap phase --reference ../../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta -o IL-SH-4_BE-OM-2.phased.bcftools.vcf IL-SH-4_BE-OM-2.bcftools.vcf ../IL-SH-4.0.08.mapped.md.bam ../BE-OM-2.0.08.mapped.md.bam
whatshap phase --reference ../../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta -o ES-DO1-1_BE-OM-2.phased.bcftools.vcf ES-DO1-1_BE-OM-2.bcftools.vcf ../ES-DO1-1.0.18.mapped.md.bam ../BE-OM-2.0.08.mapped.md.bam
whatshap phase --reference ../../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta -o RU-BAYA1-1_BE-OM-2.phased.bcftools.vcf RU-BAYA1-1_BE-OM-2.bcftools.vcf ../RU-BAYA1-1.0.07.mapped.md.bam ../BE-OM-2.0.08.mapped.md.bam
whatshap phase --reference ../../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta -o IL-SH-4_ES-DO1-1.phased.bcftools.vcf IL-SH-4_ES-DO1-1.bcftools.vcf ../IL-SH-4.0.08.mapped.md.bam ../ES-DO1-1.0.18.mapped.md.bam 
whatshap phase --reference ../../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta -o FI-BR1-39-6_RU-BAYA1-1.phased.bcftools.vcf FI-BR1-39-6_RU-BAYA1-1.bcftools.vcf ../FI-BR1-39-6.0.23.mapped.md.bam ../RU-BAYA1-1.0.07.mapped.md.bam
whatshap phase --reference ../../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta -o IL-SH-4_RU-BAYA1-1.phased.bcftools.vcf IL-SH-4_RU-BAYA1-1.bcftools.vcf ../IL-SH-4.0.08.mapped.md.bam ../RU-BAYA1-1.0.07.mapped.md.bam
whatshap phase --reference ../../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta -o ES-DO1-1_RU-BAYA1-1.phased.bcftools.vcf ES-DO1-1_RU-BAYA1-1.bcftools.vcf ../ES-DO1-1.0.18.mapped.md.bam ../RU-BAYA1-1.0.07.mapped.md.bam
whatshap phase --reference ../../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta -o BE-OM-2_BE-T1-2.phased.bcftools.vcf BE-OM-2_BE-T1-2.bcftools.vcf ../BE-OM-2.0.08.mapped.md.bam ../BE-T1-2.0.35.mapped.md.bam
whatshap phase --reference ../../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta -o BE-OM-2_BE-WH2-11.phased.bcftools.vcf BE-OM-2_BE-WH2-11.bcftools.vcf ../BE-OM-2.0.08.mapped.md.bam ../BE-WH2-11.0.34.mapped.md.bam
whatshap phase --reference ../../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta -o BE-WH2-11_BE-T1-2.phased.bcftools.vcf BE-WH2-11_BE-T1-2.bcftools.vcf ../BE-T1-2.0.35.mapped.md.bam ../BE-WH2-11.0.34.mapped.md.bam


# phasing or masking of remaining positions & generating of pseudo-diploids
for sampleduo in *.phased.bcftools.vcf
do
name="$(echo "$sampleduo" | sed -E 's/.phased.bcftools.vcf//g')"

# replace unphasable sites either with reference or alternative allele (0/1 or 1/0 -> reference; 1/1 -> alternative)
sed -e $'s~0/1:[^\t]*\t~0/0:0\t~g' $sampleduo | sed -e $'s~1/0:[^\t]*\t~0/0:0\t~g' | sed -e $'s~\t0/1:[^$]*$~\t0/0:0~g' | sed -e $'s~\t1/0:[^$]*$~\t0/0:0~g' | sed -e 's~1/1~1|1~g' | sed -e 's~0/0~0|0~g' > $name.allphased.bcftools.vcf

# combine the haplotypes
# previously heterozygous sites are being masked with N (including such, which were replaced above and not heterozygous when combined)
# 1_a 2_a
sed -E $'s~([0-1])\|[0-1]:[^\t]*\t([0-1]).*~\\1|\\2~g' $name.allphased.bcftools.vcf | awk '{if($10 == "0|0")$5="."} {if($10 == "1|1")$4="."} {print $0}' | vcfutils.pl vcf2fq -d 16 -D 100 | gzip > $name.1a2a.fq.gz
# 1_a 2_b
sed -E $'s~([0-1])\|[0-1]:[^\t]*\t[0-1]\|([0-1]).*~\\1|\\2~g' $name.allphased.bcftools.vcf | awk '{if($10 == "0|0")$5="."} {if($10 == "1|1")$4="."} {print $0}' | vcfutils.pl vcf2fq -d 16 -D 100 | gzip > $name.1a2b.fq.gz
# 1_b 2_a
sed -E $'s~[0-1]\|([0-1]):[^\t]*\t([0-1]).*~\\1|\\2~g' $name.allphased.bcftools.vcf | awk '{if($10 == "0|0")$5="."} {if($10 == "1|1")$4="."} {print $0}' | vcfutils.pl vcf2fq -d 16 -D 100 | gzip > $name.1b2a.fq.gz
# 1_b 2_b
sed -E $'s~[0-1]\|([0-1]):[^\t]*\t[0-1]\|([0-1]).*~\\1|\\2~g' $name.allphased.bcftools.vcf | awk '{if($10 == "0|0")$5="."} {if($10 == "1|1")$4="."} {print $0}' | vcfutils.pl vcf2fq -d 16 -D 100 | gzip > $name.1b2b.fq.gz

done

# psmc commands
for pseudodiploid in *.fq.gz
do
name="$(echo "$pseudodiploid" | sed -E 's/.fq.gz//g')"
echo $name
../../../../../../bioinformatics/psmc/psmc-master/utils/fq2psmcfa -q20 $name.fq.gz > $name.psmcfa
../../../../../../bioinformatics/psmc/psmc-master/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o $name.psmc $name.psmcfa
../../../../../../bioinformatics/psmc/psmc-master/utils/psmc_plot.pl $name $name.psmc
done
