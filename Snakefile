samples = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"]
# have a list of the samples to reference throughout pipeline analysis
rule all:
    input:
        "Jetton_PipelineReport.txt"
#create a rule that all other rules reference back to and have everything feed into 
rule download_hcmv_genome:
    output:
        fasta = "genomes/hcmv.fasta"
    shell:
        """
        datasets download genome accession GCF_000845245.1 --include genome
        unzip ncbi_dataset.zip
        cp ncbi_dataset/data/GCF_000845245.1/*.fna {output.fasta}
        """
#download the HCMV reference genome using ncbi datasets
#then copy the .fna file into a genomes directory
rule bowtie2_index:
    input:
        fasta = "genomes/hcmv.fasta"
    output:
        "genomes/hcmv_index.1.bt2"
    shell:
        "bowtie2-build {input.fasta} genomes/hcmv_index"
#build a bowtie2 index from the HCMV ref genome for alignment
#I used the .1.bt2 file as a proxy to show the full index has been built
rule bowtie2_map:
    input:
        index = "genomes/hcmv_index.1.bt2",
        fwd = "downloaded_reads/{sample}_1.fastq",
        rev = "downloaded_reads/{sample}_2.fastq"
    output:
        fwd = "mapped/{sample}_1.fastq",
        rev = "mapped/{sample}_2.fastq"
    shell:
        """
        before=$(( $(wc -l < {input.fwd}) / 4 ))
        bowtie2 -x genomes/hcmv_index -1 {input.fwd} -2 {input.rev} --al-conc mapped/{wildcards.sample}_%.fastq --quiet
        after=$(( $(wc -l < {output.fwd}) / 4 ))
        echo 'Sample {wildcards.sample} had $before read pairs before and $after read pairs after Bowtie2 filtering.' >> bowtie2_report.txt
        """
#aligned paired-end reads to HCMV and kept rads that successfully mapped
#read in the counts before and after filtering and log them to the report
rule spades_assembly:
    input:
        fwd = "mapped/{sample}_1.fastq",
        rev = "mapped/{sample}_2.fastq"
    output:
        contigs = "assemblies/{sample}/contigs.fasta"
    shell:
        """
        spades.py -k 99 -t 2 --only-assembler -1 {input.fwd} -2 {input.rev} -o assemblies/{wildcards.sample}
        """
#assemble the reads that mapped to the hcmv genome using spades
#only using a k-mer size of 99 and skipping the read error correction step
rule contig_stats:
    input:
        contigs = "assemblies/{sample}/contigs.fasta"
    output:
        longest = "longest_contigs/{sample}_longest.fasta"
    script:
        "scripts/contig_finder.py"
#runs a python script that extracts the # of contigs >1000 bp and overall length
#also extracts the longest kmer from each sample for use in blast
rule make_blast_db:
    output:
        blast_db = "blast_db/betaherpes.nhr"
    shell:
        """
        datasets download virus genome taxon "Betaherpesvirinae" --filename blast_db/betaherpes.zip
        unzip blast_db/betaherpes.zip -d blast_db
        cat blast_db/ncbi_dataset/data/*.fna > blast_db/betaherpes.fasta
        makeblastdb -in blast_db/betaherpes.fasta -dbtype nucl -out blast_db/betaherpes
        """
# download all bethaherpesvirninae genoems from NCBI Virus to build a local db
# so we can query the longest contigs against it
rule blast_query:
    input:
        query = "longest_contigs/{sample}_longest.fasta",
        blast_db = "blast_db/betaherpes.nhr"
    output:
        blast_output = "blast_results/{sample}_blast.txt"
    shell:
        "blastn -query {input.query} -db blast_db/betaherpes -out {output.blast_output} -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle' -max_hsps 1 -max_target_seqs 5 -num_threads 2"
# blast each samples longest contig against the local db
# output format 6 returns tabular results with certain fields
# the blast results are alos limited to top 5 target seqs, 1 HSP each
rule finalize_report:
    input:
        bowtie_results = "bowtie2_report.txt",
        blast_results = expand("blast_results/{sample}_blast.txt", sample=samples)
    output:
        report = "Jetton_PipelineReport.txt"
    run:
        with open(output.report, "w") as report:
            with open(input.bowtie_results) as bowtie:
                for line in bowtie:
                    report.write(line)
            for sample, blast_file in zip(samples, input.blast_results):
                report.write(f"\n{sample}: sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
                with open(blast_file) as blast:
                    for line in blast:
                        report.write(line)
# compile all blast results into the pipeline report
# expand makes the rule waits until all samples are finished in blast
# writes a header line followed by the certain fields from blast for each sample




