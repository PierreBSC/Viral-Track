# Manual for Viral-Track

Viral-Track is a R-based computational software based on **STAR** and **samtools** developped to detect and identify viruses from single-cell RNA-sequencing (scRNA-seq) raw data. This tool was tested on various scRNA-seq datasets including mouse and human infected samples as described in our paper *'Detecting and studying viral infection at the single-cell resolution using Viral-Track'*. 


Installation
-------------

Before running Viral-Track, several dependencies must be installed :

1 . The first step is to install [**R software**] (https://www.r-project.org/). Once this is done, several Bioconductor packages  have to be installed too. To do so start a R session and type :


```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.10") 
BiocManager::install(c("Biostrings", "ShortRead","doParallel","GenomicAlignments"))
```


2 . **S**pliced **T**ranscript **A**lignment to **A** **R**eference (**STAR**) has to be installed. The full installation process is described on the STAR [Github](https://github.com/alexdobin/STAR). On Ubuntu this can be done directly by typing :

```batch
sudo apt install rna-star
```
3 . **Samtools** suite is also required. The installation is described extensively [here](http://www.htslib.org/download/). On Ubuntu this is done by simply typing :

```batch
sudo apt install samtools
```

4 . Lastly, the transcript assembler **StringTie** is needed. This installation process is described [here](https://ccb.jhu.edu/software/stringtie/). Don't forget to add StringTie to your shell's PATH directory.

Creation of the  Index 
----------

The first step consists in creating a **STAR** index that include both host and virus reference genomes.
To do so first download the **ViruSite**  [genome reference database](http://www.virusite.org/index.php?nav=download). Host genome has also to be downloaded from the [ensembl website](https://www.ensembl.org/info/data/ftp/index.html). 

The  **STAR** index can now be build by typing :

```batch
mkdir /path/to/index
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /path/to/index --genomeFastaFiles /path/to/Virusite_file.fa  path/to/Host_genome_chromosome*.fa 
```

This can take some time and requires large amount of memory and storage space : please check that you have at least 32 GB of RAM and more than 100GB of avaible memory.

Detection of virus from scRNA-seq data
---------------

This can take some time and requires large amount of memory and storage space : please check that you have at least 32 GB of RAM and more than 100GB of avaible memory.





