# COMP383-Pipeline-Project
A project that involved assembling a Python pipeline to read in sequence data, map reads to a reference genome, assemble genomes based on mapped reads, and query the longest contigs from those assemblies against a BLAST database
Before running the pipeline, I downloaded sequence data using wget and unpacked them using fasterq. The exact script I used is here:

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045

retrieved every file using wget and the link from SRA

fasterq-dump ./SRR5660030
fasterq-dump ./SRR5660033
fasterq-dump ./SRR5660044
fasterq-dump ./SRR5660045

used faster dump to decompress files into both paired-end reads. I then moved these sequences to a subdirectory titled "downloaded_reads".

Before running the data through the pipeline, I had to enter: "touch bowtie2_report.txt" to create an empty text file for use
