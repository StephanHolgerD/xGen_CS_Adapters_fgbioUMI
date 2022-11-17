import os
import pandas as pd
from collections import Counter
DIRECTION=["1","2"]
GROUPS=set()
sample_sheet=pd.read_csv("samples.tsv", sep="\t",dtype=object)
genome = '/home/stephano/Documents/02_gitHubProjects/DrukBam_Variant/test_data/chr7.fasta'


def check_symlink(file1, file2):
    try:
        os.symlink(file1, file2)
    except FileExistsError:
        print("Link for " + file1 + " is already present in 01_raw")

SAMPLES=list(sample_sheet["ID"])


for b in set(SAMPLES):
    os.makedirs("../01_raw/" +b+ "/fastqc", exist_ok=True)
    check_symlink(sample_sheet[sample_sheet["ID"]==b]["forward reads"].values[0], "../01_raw/"+b +"/"+ b +"_1P.fastq.gz")
    check_symlink(sample_sheet[sample_sheet["ID"]==b]["reverse reads"].values[0], "../01_raw/"+b +"/"+ b +"_2P.fastq.gz")


rule all:
    input:
        expand('../01_raw/{sample}/fastqc/{sample}_{direction}P_fastqc.html',direction=DIRECTION, sample=SAMPLES),
        expand('../04_ConsensusMappedBam/{sample}.consensus.bam', sample=SAMPLES)

rule fastqc1:
    input:
        r = '../01_raw/{sample}/{sample}_{direction}P.fastq.gz',
    threads: 2
    priority: 50
    output:
        '../01_raw/{sample}/fastqc/{sample}_{direction}P_fastqc.html'
    conda:
        "envs/fastqc.yaml"
    shell:
        'fastqc -o ../01_raw/{wildcards.sample}/fastqc -t {threads} --extract {input.r}'



rule Fastq2Sam:
    input:
        r1= '../01_raw/{sample}/{sample}_1P.fastq.gz',
        r2= '../01_raw/{sample}/{sample}_2P.fastq.gz',
    output:
        unmappedBam = '../02_unmappedBam/{sample}.unmapped.bam'
        
    conda:
        "envs/picard_bwa.yaml"
        
    threads:1
    
    shell:
        'picard FastqToSam FASTQ={input.r1} FASTQ2={input.r2} O={output.unmappedBam} SM=sample'
        



rule ExtractUmisFromBam:
    input:
        unmappedBam = '../02_unmappedBam/{sample}.unmapped.bam'
    output:
        unmappedBamwithUMI = '../02_unmappedBam/{sample}.unmappedUMI.bam'
    conda:
        "envs/fgbio.yaml"
    threads:1
    shell:
        'fgbio ExtractUmisFromBam --input={input.unmappedBam} --output={output.unmappedBamwithUMI} --read-structure=3M2S146T 3M2S146T --molecular-index-tags=ZA ZB --single-tag=RX'




rule Index:
    input:
        unmappedBamwithUMI = '../02_unmappedBam/{sample}.unmappedUMI.bam'
    output:
        unmappedBamwithUMIbai = '../02_unmappedBam/{sample}.unmappedUMI.bam.bai'
    conda:
        "envs/samtools.yaml"
    threads:1
    shell:
        'samtools index {input.unmappedBamwithUMI}'
    




rule Mapping:
    input:
        unmappedBamwithUMIbai = '../02_unmappedBam/{sample}.unmappedUMI.bam.bai',
        unmappedBamwithUMI = '../02_unmappedBam/{sample}.unmappedUMI.bam'

    output:
        mappedBAM = '../03_mappedBam/{sample}.mapped.bam'

    conda:
        "envs/picard_bwa.yaml"
    threads:4
    
    shell:
        'picard SamToFastq I={input.unmappedBamwithUMI} F=/dev/stdout INTERLEAVE=true \
            | bwa mem -p -t {threads} {genome} /dev/stdin \
                | picard MergeBamAlignment \
                    UNMAPPED={input.unmappedBamwithUMI} ALIGNED=/dev/stdin O={output.mappedBAM} R={genome} \
                    SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 \
                    ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true'
                    
rule GroupReadsByUmi:
    input:
        mappedBAM = '../03_mappedBam/{sample}.mapped.bam'
    output:
        groupedmappedBAM = '../03_mappedBam/{sample}.grouped.mapped.bam'
     
    threads:1
       
    conda:
        "envs/fgbio.yaml"

    shell:
        'fgbio GroupReadsByUmi --input={input.mappedBAM} --output={output.groupedmappedBAM} --strategy=paired --edits=1 --min-map-q=2'




rule CallDuplexConsensusReads:
    input:
        groupedmappedBAM = '../03_mappedBam/{sample}.grouped.mapped.bam'
        
    output:
        groupedmappedconsensusBAM = '../03_mappedBam/{sample}.grouped.mapped.consensus.bam'

    threads:1
    
    conda:
        "envs/fgbio.yaml"

    shell:
        'fgbio CallDuplexConsensusReads \
            --input={input.groupedmappedBAM} --output={output.groupedmappedconsensusBAM} \
            --error-rate-pre-umi=45 --error-rate-post-umi=30 \
            --min-input-base-quality=30'


rule CallMappConsensusReads:
    input:
        groupedmappedconsensusBAM = '../03_mappedBam/{sample}.grouped.mapped.consensus.bam'
        
    output:
        consensusmapped2BAM= '../04_ConsensusMappedBam/{sample}.consensus.bam'
        
    threads:4

    conda:
        "envs/picard_bwa.yaml"

    shell:
        'picard SamToFastq I={input.groupedmappedconsensusBAM} \
            F=/dev/stdout INTERLEAVE=true \
            | bwa mem -p -t {threads} {genome} /dev/stdin \
            | picard MergeBamAlignment \
            UNMAPPED={input.groupedmappedconsensusBAM}  ALIGNED=/dev/stdin \
            O={output.consensusmapped2BAM} R={genome} \
            SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 \
            ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true'