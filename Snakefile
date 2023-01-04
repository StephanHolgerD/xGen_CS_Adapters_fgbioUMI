import os
import pandas as pd
from collections import Counter

#change to the path of your ref genome
genome = '/mnt/d/2022/MOR/MOR_22_reanalyseTechnicalDataset/00_ref/chr7.fasta'
tmp = '/home/drukewitz/conda_tmp/'


bed = 'data/Twist-custom_hg19.bed'






DIRECTION=["1","2"]
GROUPS=set()
sample_sheet=pd.read_csv("samples.tsv", sep="\t",dtype=object)


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
        expand('../06_Vardict/{sample}.vcf', sample=SAMPLES),




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
        'picard FastqToSam TMP_DIR={tmp} FASTQ={input.r1} FASTQ2={input.r2} O={output.unmappedBam} SM=sample'
        



rule ExtractUmisFromBam:
    input:
        unmappedBam = '../02_unmappedBam/{sample}.unmapped.bam'
    output:
        unmappedBamwithUMI = '../02_unmappedBam/{sample}.unmappedUMI.bam'
    conda:
        "envs/fgbio.yaml"
    threads:1
    shell:
        'fgbio -XX:-UseGCOverheadLimit -Xms750m -Xmx24g --tmp-dir={tmp} ExtractUmisFromBam --input={input.unmappedBam} --output={output.unmappedBamwithUMI} --read-structure=3M2S146T 3M2S146T --molecular-index-tags=ZA ZB --single-tag=RX'



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
    threads:12
    
    shell:
        'picard SamToFastq TMP_DIR={tmp} I={input.unmappedBamwithUMI} F=/dev/stdout INTERLEAVE=true \
            | bwa mem -p -t {threads} {genome} /dev/stdin \
                | picard MergeBamAlignment  \
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
        'fgbio -XX:-UseGCOverheadLimit -Xms750m -Xmx24g --tmp-dir={tmp} GroupReadsByUmi --input={input.mappedBAM} --output={output.groupedmappedBAM} --strategy=paired --edits=1 --min-map-q=2'



###################################################################################################################################################


rule CallMolecularConsensusReads:
    input:
        groupedmappedBAM = '../03_mappedBam/{sample}.grouped.mapped.bam'
        
    output:
        groupedUnmappedconsensusBAM = '../04_CallMolecularConsensusReads/{sample}.grouped.mapped.consensus.bam'

    threads:1
    
    conda:
        "envs/fgbio.yaml"

    shell:
        'fgbio -XX:-UseGCOverheadLimit \
            -Xms750m -Xmx24g --tmp-dir={tmp} \
            CallMolecularConsensusReads --min-input-base-quality=10 --min-reads=2 --max-reads=1000000 --output-per-base-tags=false --sort-order=:none: -i {input.groupedmappedBAM} -o {output.groupedUnmappedconsensusBAM}'


###################################################################################################################################################



rule SortBam:
    input:
        groupedUnmappedconsensusBAM = '../04_CallMolecularConsensusReads/{sample}.grouped.mapped.consensus.bam'
    output:
        groupedUnmappedconsensusBAMsorted = '../04_CallMolecularConsensusReads/{sample}.grouped.mapped.consensus.sorted.bam'
    threads:1
    
    conda:
        "envs/picard_bwa.yaml"

    shell:
        "picard SortSam -I {input.groupedUnmappedconsensusBAM} -O {output.groupedUnmappedconsensusBAMsorted} -SO queryname"




rule CallMappConsensusReads:
    input:
        groupedmappedconsensusBAM = '../04_CallMolecularConsensusReads/{sample}.grouped.mapped.consensus.sorted.bam'
        
    output:
        consensusmappedBAM= '../05_ConsensusMappedBam/{sample}.consensus.bam'
        
    threads:12

    conda:
        "envs/picard_bwa.yaml"

    shell:
        'picard SamToFastq TMP_DIR={tmp} I={input.groupedmappedconsensusBAM} \
            F=/dev/stdout INTERLEAVE=true \
            | bwa mem -p -t {threads} {genome} /dev/stdin \
            | picard MergeBamAlignment \
            UNMAPPED={input.groupedmappedconsensusBAM}  ALIGNED=/dev/stdin O={output.consensusmappedBAM} R={genome} \
            SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 \
            ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true'


rule VardictCalling:
    input:
        consensusmappedBAM= '../05_ConsensusMappedBam/{sample}.consensus.bam'
    
    threads: 6

    output:
        vcf= '../06_Vardict/{sample}.vcf'

    conda:
        "envs/vardict.yaml"

    shell:
        'vardict-java -G {genome} \
            -th {threads} -f 0.001 -r 2 \
            -b {input.consensusmappedBAM} \
            -z -c 1 -S 2 -E 3 -g 4 \
            --nosv {bed}| teststrandbias.R | var2vcf_valid.pl -E -f 0.001 > {output.vcf}'