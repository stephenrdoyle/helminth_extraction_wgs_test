# DNA extraction and whole genome sequencing of single egg and larval samples from helminths

## Table of contents
0. [Project overview](#overview)
1. [Projects setup](#setup)
2. [Mapping](#mapping)


## 00 - Project overview

The aim of this work is to....





## 01 - Project setup <a name="setup"></a>
### Setup a working environment for the analysis.
``` shell
mkdir $HOME/HELMINTH_EXTRACTION_WGS
cd $HOME/HELMINTH_EXTRACTION_WGS

# make working directories
mkdir 00_SCRIPTS 01_REFERENCES 02_RAWDATA 03_MAPPING 04_ANALYSIS
```

### Reference genomes
Get the reference genomes for mapping
```
cd 01_REFERENCES

# a_canium
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS12/species/ancylostoma_caninum/PRJNA72585/ancylostoma_caninum.PRJNA72585.WBPS12.genomic.fa.gz

# d_medinensis
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS12/species/dracunculus_medinensis/PRJEB500/dracunculus_medinensis.PRJEB500.WBPS12.genomic.fa.gz

# h_contortus
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS12/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS12.genomic.fa.gz

# s_mansoni
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS12/species/schistosoma_mansoni/PRJEA36577/schistosoma_mansoni.PRJEA36577.WBPS12.genomic.fa.gz

# a_dissimilis - note, we dont have a reference genome for this species, so using ascaris_lumbricoides
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS12/species/ascaris_lumbricoides/PRJEB4950/ascaris_lumbricoides.PRJEB4950.WBPS12.genomic.fa.gz

# d_immitis
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS12/species/dirofilaria_immitis/PRJEB1797/dirofilaria_immitis.PRJEB1797.WBPS12.genomic.fa.gz

# s_stercoralis
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS12/species/strongyloides_stercoralis/PRJEB528/strongyloides_stercoralis.PRJEB528.WBPS12.genomic.fa.gz

# t_muris
wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS12/species/trichuris_muris/PRJEB126/trichuris_muris.PRJEB126.WBPS12.genomic.fa.gz

# unzip
for i in *.gz; do gunzip ${i} ;  done

```

### Sample and lane list
I use this format for a few different steps in my workflows. It is simple a tab delimited file containing the sample name, and the sequencing lane ID. Samples with more
than a single sequencing lane will be represented twice, but overall, each line of the file will be unique. At Sanger, we have a pathfind script that can be looped through to retrieve the IDs, and then these can be renames as appropriate.

Example is:
Hc_L1_chip_1	21094_1#1
Hc_L1_chip_2	21094_1#2
Hc_L1_chip_3	21094_1#3
Hc_Egg_chip_1	21094_1#4
Hc_Egg_chip_2	21094_1#5
Hc_Egg_chip_3	21094_1#6


### Raw data
Get the raw sequencing data for the analysis. All data will be available for download from ENA, links to which are provided in the metadata table in which the samples are described.
Be aware the the steps in the next section may not work exactly as stated, as it is based on how I retrived the Sanger environment. However,
at the end of this subsection, it should be pretty clear how the data needs to be formatted for the mapping steps.

```shell
