# DNA extraction and whole genome sequencing of single egg and larval samples from helminths

## Table of contents
0. [Project overview](#overview)
1. [Project setup](#setup)
2. [Mapping](#mapping)
3. [Analysis](#analysis)


## 00 - Project overview

The aim of this work is to....





## 01 - Project setup <a name="setup"></a>
### Setup a working environment for the analysis.
In my work environment, it is: ${WORKING_DIR}

``` shell

mkdir ${HOME}/HELMINTH_EXTRACTION_WGS
cd ${HOME}/HELMINTH_EXTRACTION_WGS
WORKING_DIR=${PWD}

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
than a single sequencing lane will be represented twice, but overall, each line of the file will be unique. At Sanger, we have a "pathfind" script that can be looped through to retrieve the IDs, and then these can be renames as appropriate.

```shell
#Example is:
DM_LS1_NXT_LCM_001      26924_3#1
DM_LS1_NXT_LCM_002      26924_3#2
DM_LS1_NXT_LCM_003      26924_3#3
DM_LS1_NXT_LCM_004      26924_3#4
DM_LS1_NXT_LCM_005      26924_3#5
DM_LS1_NXT_LCM_006      26924_3#6
DM_LS1_NXT_LCM_007      26924_3#7
DM_LS1_NXT_LCM_008      26924_3#8
DM_LS1_NXT_LCM_009      26924_3#9
DM_LS1_NXT_LCM_010      26924_3#10
```

The actual "samples_lanes.list" used in this study is located here: XXXX

${WORKING_DIR}/02_RAWDATA/samples_lanes.list

Need to generate a per species sample_lane list for mapping, simply becasue they need to be mapped to different reference genomes. Otherwise, would have used the orginal file.
```shell
# use grep to extract lines associated with each species
for i in DM AC HC SM AD DI SS TM; do
          grep ^${i} samples_lanes.list > ${i}_samples_lanes.list;
     done
```



### Raw data
Get the raw sequencing data for the analysis. All data will be available for download from ENA, links to which are provided in the metadata table in which the samples are described.
Be aware the the steps in the next section may not work exactly as stated, as it is based on how I retrieved the Sanger environment. However,
at the end of this subsection, it should be pretty clear how the data needs to be formatted for the mapping steps.

```shell
# getting data form iRODs
kinit  # this requires a password

icd /seq/28536
ils | grep "28536.*.cram$" | grep -v "phix" | while read -r data; do iget /seq/28536/$data ./; done &

icd /seq/25934
ils | grep "25934.*.cram$" | grep -v "phix" | while read -r data; do iget /seq/25934/$data ./; done &


icd /seq/26790
ils | grep "26790.*.cram$" | grep -v "phix" | while read -r data; do iget /seq/26790/$data ./; done &

# get the old "chip" and "nextera" data - it is in pathfind
for i in 21094_1#1 21094_1#2 21094_1#3 21094_1#4 21094_1#5 21094_1#6 21033_1#1 21033_1#2 21033_1#3 21033_1#4 21033_1#5 21033_1#6 ; do pathfind -t lane -i ${i} --symlink ./ --filetype fastq ; done &

# get the GW data
while read NAME; do pathfind -t lane -i ${NAME} --symlink ./ --filetype fastq ; done < gw.list

# get the SM data
while read NAME LANE; do pathfind -t lane -i ${LANE} --symlink ./ --filetype fastq ; done < sm.list

#--- convert crams to fastq, and then gzip
for i in *cram; do samtools view -ub --threads 4 ${i} | samtools sort -n - | samtools fastq -1 ${i%.cram}_1.fastq.gz -2 ${i%.cram}_2.fastq.gz - ; done

rm *.cram

```




## 02 - Mapping to the genome <a name="mapping"></a>
Mapping uses a wrapper script called "run_bwamem_splitter" to map reads using bwa mem to the genome. It speeds up the process by splitting reads into smaller chunks, maps them individually, and then merges the mapping files.

It is quite specific to the Sanger environment - ie, samtools, bwa, and that it automatically submits jobs to out LSF system - and probably wont work anywhere else.

```shell
cd ${WORKING_DIR}/03_MAPPING

screen

while read NAME LANE; do \
~sd21/bash_scripts/run_bwamem_splitter ${NAME}_${LANE} \
${WORKING_DIR}/01_REFERENCES/dracunculus_medinensis.PRJEB500.WBPS12.genomic.fa \
${WORKING_DIR}/02_RAWDATA/${LANE}_1.fastq.gz \
${WORKING_DIR}/02_RAWDATA/${LANE}_2.fastq.gz; \
done < ${WORKING_DIR}/02_RAWDATA/DM_samples_lanes.list &

sleep 5

while read NAME LANE; do \
~sd21/bash_scripts/run_bwamem_splitter ${NAME}_${LANE} \
${WORKING_DIR}/01_REFERENCES/ancylostoma_caninum.PRJNA72585.WBPS12.genomic.fa \
${WORKING_DIR}/02_RAWDATA/${LANE}_1.fastq.gz \
${WORKING_DIR}/02_RAWDATA/${LANE}_2.fastq.gz; \
done < ${WORKING_DIR}/02_RAWDATA/AC_samples_lanes.list &

sleep 5

while read NAME LANE; do \
~sd21/bash_scripts/run_bwamem_splitter ${NAME}_${LANE} \
${WORKING_DIR}/01_REFERENCES/haemonchus_contortus.PRJEB506.WBPS12.genomic.fa \
${WORKING_DIR}/02_RAWDATA/${LANE}_1.fastq.gz \
${WORKING_DIR}/02_RAWDATA/${LANE}_2.fastq.gz; \
done < ${WORKING_DIR}/02_RAWDATA/HC_samples_lanes.list &

sleep 5

while read NAME LANE; do \
~sd21/bash_scripts/run_bwamem_splitter ${NAME}_${LANE} \
${WORKING_DIR}/01_REFERENCES/schistosoma_mansoni.PRJEA36577.WBPS12.genomic.fa \
${WORKING_DIR}/02_RAWDATA/${LANE}_1.fastq.gz \
${WORKING_DIR}/02_RAWDATA/${LANE}_2.fastq.gz; \
done < ${WORKING_DIR}/02_RAWDATA/SM_samples_lanes.list &

sleep 5

while read NAME LANE; do \
~sd21/bash_scripts/run_bwamem_splitter ${NAME}_${LANE} \
${WORKING_DIR}/01_REFERENCES/ascaris_lumbricoides.PRJEB4950.WBPS12.genomic.fa \
${WORKING_DIR}/02_RAWDATA/${LANE}_1.fastq.gz \
${WORKING_DIR}/02_RAWDATA/${LANE}_2.fastq.gz; \
done < ${WORKING_DIR}/02_RAWDATA/AD_samples_lanes.list &

sleep 5

while read NAME LANE; do \
~sd21/bash_scripts/run_bwamem_splitter ${NAME}_${LANE} \
${WORKING_DIR}/01_REFERENCES/dirofilaria_immitis.PRJEB1797.WBPS12.genomic.fa \
${WORKING_DIR}/02_RAWDATA/${LANE}_1.fastq.gz \
${WORKING_DIR}/02_RAWDATA/${LANE}_2.fastq.gz; \
done < ${WORKING_DIR}/02_RAWDATA/DI_samples_lanes.list &

sleep 5

while read NAME LANE; do \
~sd21/bash_scripts/run_bwamem_splitter ${NAME}_${LANE} \
${WORKING_DIR}/01_REFERENCES/strongyloides_stercoralis.PRJEB528.WBPS12.genomic.fa \
${WORKING_DIR}/02_RAWDATA/${LANE}_1.fastq.gz \
${WORKING_DIR}/02_RAWDATA/${LANE}_2.fastq.gz; \
done < ${WORKING_DIR}/02_RAWDATA/SS_samples_lanes.list &

sleep 5

while read NAME LANE; do \
~sd21/bash_scripts/run_bwamem_splitter ${NAME}_${LANE} \
${WORKING_DIR}/01_REFERENCES/trichuris_muris.PRJEB126.WBPS12.genomic.fa \
${WORKING_DIR}/02_RAWDATA/${LANE}_1.fastq.gz \
${WORKING_DIR}/02_RAWDATA/${LANE}_2.fastq.gz; \
done < ${WORKING_DIR}/02_RAWDATA/TM_samples_lanes.list &

# check to see there are flagstat data per mapping. If not, print name
for i in *out; do
          if [ ! -f "${i}/${i%_bwasplitter_out}.merged.sorted.marked.flagstat" ]; then
                    echo ${i};
               fi;
          done

```







## 03 - Analysis <a name="analysis"></a>

The run_bwamem_splitter script automatically determiend s flagstat and bamstats for the bam files in each mapping job. We'll summarise all of that informaiton using multiqc

```shell
cd ${WORKING_DIR}/04_analysis

multiqc ${WORKING_DIR}/03_MAPPING
