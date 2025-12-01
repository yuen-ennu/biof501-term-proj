# Nextflow Pipeline for Viral Taxonomic Classification for Short-Read Whole Genome Sequencing in Human Metagenomics.
### By: Yasmine Lau

***

## Background and Rationale

The human microbiome is the collective genome of microorganisms, including bacteria, archaea, fungi, eukaryotes, and viruses. In recent years, it has become clear that the microbiome profoundly influences human homeostasis, health and disease states, metabolism, immune function, and even neurobehavioral traits. Our gut microbiome alone contains trillions of bacteria representing thousands of species, possibly beginning to develop as early as in the womb. [[1]](#1)[[2]](#2)[[3]](#3)

Bacterial contributions to our microbiome, or bacteriome, and their effects on homeostasis have been extensively studied. Yet, the viral component (virome) has been more challenging to characterize, largely due to difficulties in detection and discovery, including the lack of viral culture methods. However, ongoing technological and metagenomic advances have begun to illuminate the role of the virome in human biology. [[4]](#4)[[5]](#5) Viral integration into the human genome can have profound effects on health and disease, impacting homeostasis and contributing to up to 10–15% of various cancers, including lymphomas, herpes-associated malignancies, cervical cancer, and other virus-associated tumors.[[6]](#6)

Whole genome shotgun sequencing has become a staple in metagonomics, able to fragment collective genomic material of all organisms in a sample, and employing bioinformatics approaches, map them to reference genomes. Alignment to the correct genomes across biological kingdoms, is certainly a key step in informative taxonomic analyses. 

We present VIRUSmap, a nextflow pipeline that identifies viral taxonomy from human metagenomic samples. The VIRUSmap pipeline is an easy to use, non-memory intensive collection of versatile and fast bioinformatics tools, made for raw, paired-end, short-read whole genome sequencing data.  
We start with raw metagenomic sequence data from NCBI’s Sequence Read Archive [(SRA)](https://www.ncbi.nlm.nih.gov/sra) database and directly input the SRA ID. This will input all sample runs into the pipeline from the same BioSample ID. If there are runs that are not paired-end, they will be excluded. 
We first use [(Fastp)](https://github.com/OpenGene/fastp) for read pre-processing, then [(fastqc)](https://github.com/s-andrews/FastQC) to perform quality control. We then align these reads with the versatile aligner [(Minimap2)](https://github.com/lh3/minimap2) ideal for whole-genome sequencing, and for our viral genomes are sourced by NCBI’s [(RefSeq viral genome collection)](https://benlangmead.github.io/aws-indexes/k2). Moving forward, we only take the unmapped reads, assuming these reads are enriched for viral sequences compared to the original metagenomic data, which we perform the viral taxonification and classification on. To do so, we use [(Kraken2)](https://github.com/DerrickWood/kraken2) and [(Bracken)](https://github.com/jenniferlu717/Bracken) for rapid taxonomic classification and abundance estimation respectively, allowing us to identify the viral sequences. With [(KrakenTools)](https://github.com/jenniferlu717/KrakenTools) we enable downstream processing for visualization with [(Krona)](https://github.com/marbl/Krona/wiki).


## Example 
VIRUSmap.nf is preloaded with an example from the [(American Gut Project)](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB11419). BioSample:[(SAMEA8947847)](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMEA8947847&o=bytes_l%3Aa). 

Below are the runs performed on this sample. The AMPLICON Assay Type is automatically filtered from the VIRUSmap pipeline.
visuals/multiple_assay_handling_example.png





## VIRUSmap dependencies
To use VIRUSmap, the user will need to have [(conda)](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and [(git)](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) installed.

## Installtaion
1. Deactivate any existing conda environemnt. 
```
conda deactivate
```
2. Clone VIRUSmap repository.
```
git clone https://github.com/yuen-ennu/biof501-term-proj
```
3. Go to project directory.
```
cd biof501-term-proj
```
4. Create VIRUSmap conda environment and activate it.
```
conda env create -f VIRUSmap.yml
conda activate VIRUSMAP
```
5. Run nextflow file.
```
nextflow run VIRUSmap.nf -c nextflow.config
```

## Outputs (see )




## References
<a id="1">[1]</a>
Valdes AM, Walter J, Segal E, Spector TD. Role of the gut microbiota in nutrition and health. BMJ. 2018;361(361):k2179.
‌
<a id="2">[2]</a>
Walker RW, Clemente JC, Peter I, Loos RJF. The prenatal gut microbiome: are we colonized with bacteria in utero? Pediatric Obesity. 2017;12(S1):3-17.

<a id="3">[3]</a>
Gregory AC, Zablocki O, Howell A, Bolduc B, Sullivan MB. The human gut virome database. bioRxiv (Cold Spring Harbor Laboratory). Published online May 31, 2019.

<a id="4">[4]</a>
Wang D. 5 challenges in understanding the role of the virome in health and disease. Spindler KR, ed. PLOS Pathogens. 2020;16(3):e1008318.

<a id="5">[5]</a>
Stephens Z, O’Brien D, Mrunal Dehankar, Roberts LR, Iyer RK, Kocher JP. Exogene: A performant workflow for detecting viral integrations from paired-end next-generation sequencing data. PLoS ONE. 2021;16(9):e0250915-e0250915.

<a id="6">[6]</a>
Bai GH, Lin SC, Hsu YH, Chen SY. The Human Virome: Viral Metagenomics, Relations with Human Diseases, and Therapeutic Applications. Viruses. 2022;14(2):278.


nextflow run </home/ylau/scratch/biof501-term-proj/VIRUSmap.nf> -preview -with-dag VIRUSmap_dag.dot
nextflow run <VIRUSmap.nf> -with-dag VIRUSmap.dot
