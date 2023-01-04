## Snakemake for xGen_CS_Adapters_fgbioUMI chemistry

* needs a working conda installation and conda in path

* Snakemake needs to be installed in the working enviroment

* sample_working.tsv provides a example of the needed sampleseheet

* use the module like this:

```
$ cd xGen_CS_Adapters_fgbioUMI
```

the current working dirtectory needs to be the xGen_CS_Adapters_fgbioUMI folder 

```
snakemake --use-conda --cores 8 --conda-prefix "folder where conda should store the envs data (is optional)"

```


* the pipeline creates a folder structure and writes all the output in the parent directory of the current working directory 

#Folder Structure
```
.
├── 01_raw
├── 02_unmappedBam
├── 03_mappedBam
├── 04_CallMolecularConsensusReads
├── 05_ConsensusMappedBam
├── 06_Vardict
└── xGen_CS_Adapters_fgbioUMI
    ├── data
    │   ├── Twist-custom_hg19.bed
    ├── envs
    │   ├── bedtools.yaml
    │   ├── fastqc.yaml
    │   ├── mapping_viz.yaml
    │   ├── sambamba.yaml
    │   ├── samtools.yaml
    │   ├── samtools.yaml.save
    │   ├── star.yaml
    │   ├── subread.yaml
    │   └── trimmomatic.yaml
    ├── samples.tsv
    ├── sample_working.tsv
    ├── scripts
    │   ├── mapping.py
    │   └── sub_sampleVIZ.py
    └── Snakefile
```
# Author
Stephan drukewitz --> Stephan.Drukewitz@uniklinikum-dresden.de

