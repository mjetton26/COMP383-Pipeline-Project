from Bio import SeqIO

infile = snakemake.input.contigs
outfile = snakemake.output.longest
sample = snakemake.wildcards.sample


def get_contigs(input_file):
    contigs = [] #import fasta files and read them into a list
    for seq_record in SeqIO.parse(input_file, "fasta"):
        contigs.append(seq_record)

    return contigs

def get_long_contigs(contigs, min_length=1000): #return all contigs with a length longer than a 1000 bp
    return [c for c in contigs if len(c.seq) >= min_length]

def get_total_bp(contigs): #get a total length of all contigs 
    return sum(len(c.seq) for c in contigs)

all_contigs  = get_contigs(infile)
long_contigs = get_long_contigs(all_contigs, 1000)
total_bp     = get_total_bp(long_contigs)
count        = len(long_contigs)
#write to pipelinereport the statistics needed
with open("PipelineReport.txt", "a") as report:
    report.write(
        f"In the assembly of sample {sample}, "
        f"there are {count} contigs > 1000 bp and {total_bp} total bp.\n"
    )
#write the longest contigs from each sequence to a file for blast 
longest = max(all_contigs, key=lambda c: len(c.seq))
SeqIO.write(longest, outfile, "fasta")






    

