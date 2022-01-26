####### Convert Sequel generated bam file to fasta sequence needed as input for CANU assembler
### use conda environment with various pbbioconda tools installed
### only required tool here is bam2fastx (https://github.com/pacificbiosciences/bam2fastx/)

bam2fasta -o m54273_200909_122406 m54273_200909_122406.subreads.bam

### resultant file is called m54273_200909_122406.fasta.gz

####### Run CANU assembler using default rather than metagenome settings
####### Metagenome settings will come later
### use conda environment with CANU (https://canu.readthedocs.io/en/latest/) 
### and dependencies installed.

canu -p canu -d canu genomeSize=200m maxThreads=55 -pacbio-raw m54273_200909_122406.fasta.gz

### resultant assembly renamed to FI-OER-3-3.mixed_assembly.v0.1.fasta.gz
### proceed with blob analysis of genome

