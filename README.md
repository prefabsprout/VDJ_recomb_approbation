# Approbation of V(D)J Recombination Probability Model

## Goals of The Project

- To approbe V(D)J recombination probability model, which was made in previous year.
- To use results of approbation for new model.

## Methods

In this project we used a raw reads in FASTQ format (NCBI BioProject accession - PRJNA396773). Data were obtained from longitude study. Study was focused on monitoring of immune repertoires of HIV infected patient. Pair-end reads was assembled with Pandaseq tool. Assembled reads was preproccesed with Biopython - we complementary reversed reads with usage of that library. After that with Partis we cluster sequences into clonal families. Result yaml files parsed with Python script into FASTA files and with IMGT HighVQuest we annotate our immune repertoires. After that we approbe our probability model with data we obtained. With R script we analyzed results.


## How to work with our pipeline

![](/home/stephen/Git_Repositories/VDJ_Recombination/pipeline.jpg)

- To assemble raw FASTQ reads use that line: 
`pandaseq -f SRR5888724_1.fastq -r SRR5888724_2.fastq > SRR5888724.fasta`

- Data preproccesing:
To revcomp data use revcomp.py. In "default" variable insert path to your original data, in "result" path to your revcomped file. 
Also we erased corrupted string from SRR5888726 with erase_by_id.py.

- To work with Partis use next oneliners:
`sudo docker cp SRR5888729_revcomp.fasta container-1:partis/SRR5888729_revcomp.fasta # Copy FASTA file to container`
`sudo docker run -it --name container-1 -v ~:/home/stephen psathyrella/partis /bin/bash # Run Docker container with Partis`
`./bin/partis partition --infname ./SRR5888729_revcomp.fasta  --outfname ./SRR5888729_parted.yaml --n-random-queries 10000 --n-procs 2 # Run partition with sample `
`bin/partis view-output --outfname SRR5888729_parted.yaml > toparse_SRR5888729_parted.yaml # To print the partitions and/or annotations from an existing output file`
`sudo docker cp container-1:partis/toparse_SRR5888729_parted.yaml toparse_SRR5888729_parted.yaml # Copy FASTA out of container`

- To prepare our FASTA file for yaml parser use:
`sed -i 's/:/c/g' SRR5888729_revcomp.fasta`

- To parse yaml:
`python3 yaml_parser.py -i SRR5888729_revcomp.fasta -y toparse_SRR5888729_parted.yaml -o SRR5888729_parted.fasta`

- To annotate result files use http://www.imgt.org/HighV-QUEST/

- To use model put IMGT files into Model/data and launch IG.py

- To analyze results in res_analysis.r:
`experim <- read.csv('Model/experimental_data.csv') # Line 99`
`mod <- read.csv('Model/model_data.csv') # Line 100`

## Approbation results

![Гистограммы распределений длин N1-зон. Синий цвет - модельные данные, красный цвет - экспериментальные](/home/stephen/Git_Repositories/VDJ_Recombination/Model/pheno13.png)

![Таблица с результатами](/home/stephen/Git_Repositories/VDJ_Recombination/Visualisation/results_table.jpg)

## Minimum system requirements

- Ubuntu 18.04 
- Python 3.6.8
- R 3.6.1

## Materials

- https://github.com/PazhenkovaEA/VH-replacement-analysis/tree/master/VanderHeiden2017_rep_analisys
- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA396773
- Kenneth Murphy, Casey Weaver, Janeway's immunobiology, 9th edition. | New York, NY : Garland
Science/Taylor & Francis Group, LLC, 2016
- Landais E, Murrell B, Briney B, Murrell S, Rantalainen K, Berndsen Z, et al. . HIV envelope glycoform heterogeneity and localized diversity govern the initiation and maturation of a V2 apex broadly neutralizing antibody lineage. Immunity. (2017) 47:990–1003.e9. 10.1016/j.immuni.2017.11.002 
