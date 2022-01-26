####### In order to do the blobology analysis one needs to first generate taxon matches for the full contig set
### use conda environment with blast and diamond installed
### database for blast is ncbi nt; database for diamond is uniprot refseq

blastn -db /home/peter/bioinformatics/ncbi_db/nt \
       -query FI-OER-3-3.mixed_assembly.v0.1.fasta \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -out blast.out

diamond blastx \
        --query FI-OER-3-3.mixed_assembly.v0.1.fasta  \
        --db /home/peter/bioinformatics/uniprot/reference_proteomes.dmnd \
        --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-25 \
        > diamond.out

### we can also take a look at the coverage of individual contigs from the PacBio data used for the assembly
### need a conda environment with minimap2 and samtools

minimap2 -ax map-pb -t 55 FI-OER-3-3.mixed_assembly.v0.1.fasta m54273_200909_122406.fasta.gz | \
samtools sort -O BAM -o FI-OER-3-3.align.bam -


### now we can load these hits and coverage into a blobtools project
### need a conda environment with blobtools v.1.0+ installed

blobtools create -i I-OER-3-3.mixed_assembly.v0.1.fasta -b FI-OER-3-3.align.bam \
 -t blast.out -t diamond.out -o FI-OER-3-3.blobtools
blobtools view -i FI-OER-3-3.blobtools.blobDB.json
blobtools plot -i FI-OER-3-3.blobtools.blobDB.json

