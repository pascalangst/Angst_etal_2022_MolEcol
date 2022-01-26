####### This PSMC analysis is used to reconstruct the demographic history 

### The first step is a downsampling of all samples to a common average coverage of 50X

samtools view -s 123.08 BE-OM-2.mapped.md.bam -o ../../psmc/BE-OM-2.0.08.mapped.md.bam
samtools view -s 123.23 FI-BR1-39-6.mapped.md.bam -o ../../psmc/FI-BR1-39-6.0.23.mapped.md.bam
samtools view -s 123.08 IL-SH-4.mapped.md.bam -o ../../psmc/IL-SH-4.0.08.mapped.md.bam
samtools view -s 123.35 BE-T1-2.mapped.md.bam -o ../../psmc/BE-T1-2.0.35.mapped.md.bam
samtools view -s 123.34 BE-WH2-11.mapped.md.bam -o ../../psmc/BE-WH2-11.0.34.mapped.md.bam
samtools view -s 123.08 CY-PA3-2.mapped.md.bam -o ../../psmc/CY-PA3-2.0.08.mapped.md.bam
samtools view -s 123.18 ES-DO1-1.mapped.md.bam -o ../../psmc/ES-DO1-1.0.18.mapped.md.bam
samtools view -s 123.14 ES-OM-5.mapped.md.bam -o ../../psmc/ES-OM-5.0.14.mapped.md.bam
samtools view -s 123.53 FI-FHS2-11-3.mapped.md.bam -o ../../psmc/FI-FHS2-11-3.0.53.mapped.md.bam
samtools view -s 123.05 FI-OER-3-3.mapped.md.bam -o ../../psmc/FI-OER-3-3.0.05.mapped.md.bam
samtools view -s 123.07 IL-G-3.mapped.md.bam -o ../../psmc/IL-G-3.0.07.mapped.md.bam
samtools view -s 123.07 IL-PS-6.mapped.md.bam -o ../../psmc/IL-PS-6.0.07.mapped.md.bam
samtools view -s 123.07 RU-BAYA1-1.mapped.md.bam -o ../../psmc/RU-BAYA1-1.0.07.mapped.md.bam
samtools view -s 123.13 RU-TY1-3.mapped.md.bam -o ../../psmc/RU-TY1-3.0.13.mapped.md.bam
samtools view -s 123.10 RU-TY2-2.mapped.md.bam -o ../../psmc/RU-TY2-2.0.10.mapped.md.bam
samtools view -s 123.13 RU-ZB1-1.mapped.md.bam -o ../../psmc/RU-ZB1-1.0.13.mapped.md.bam
samtools view -s 123.09 SE-G4-20.mapped.md.bam -o ../../psmc/SE-G4-20.0.09.mapped.md.bam
samtools view -s 123.15 SE-H1-4.mapped.md.bam -o ../../psmc/SE-H1-4.0.15.mapped.md.bam

cd ../../psmc/


### This is the actual code for running PSMC

for downsampled in *bam
do
name="$(echo "$downsampled" | sed -E 's/.mapped.*//g')"
echo $name
samtools mpileup -C50 -uf ../FI-OER-3-3.pacbio.v1/FI-OER-3-3.pacbio.v1.fasta $downsampled | bcftools call -c - \ | vcfutils.pl vcf2fq -d 16 -D 100 | gzip > $name.fq.gz
./psmc/psmc-master/utils/fq2psmcfa -q20 $name.fq.gz > $name.psmcfa
./psmc/psmc-master/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o $name.psmc $name.psmcfa
.psmc/psmc-master/utils/psmc_plot.pl $name $name.psmc
done

### This will output one .psmc file per sample. The plotting is done with the according R script
