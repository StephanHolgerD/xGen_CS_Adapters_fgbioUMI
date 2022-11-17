from Bio import SeqIO
import sys
import gzip
a=0
with gzip.open(sys.argv[1], "rt") as fastq:
    for record in SeqIO.parse(fastq, "fastq"):
        a=a+1
print(a)
