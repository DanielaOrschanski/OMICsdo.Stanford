# OMICsdo.Stanford

R package for RNAseq data analysis.

Using the command **OMICSdo.Stanford::setup_environment()** it downloads the following tools:
1. FastQC
2. Trim Galore
3. STAR
4. Arriba
5. R packages required:
   GEOquery, Rsubread, ggplot2, ggpubr, openxlsx, png, readr, readxl, reshape2, rmarkdown, rstudioapi, stringr, survival, survminer, viridis, qckitfastq.

Additionally, it downloads the hg38 genome reference (fasta and annotation).

**Important** As soon as setup_environment() is executed, a window will pop up and you will have to indicate the folder that will store all those resources. At lest 30 GB is required.

System installations are needed:

sudo apt install wget unzip cutadapt tar bzip2

## Practical Steps:

The fastq files should be organized in the following structure:

/Samples

--  /Samples/ID

---     /Samples/ID/ID_R1.fastq

---     /Samples/ID/ID_R2.fastq
            

### 1. Instalation of OMICsdo.Stanford from Github:

library(devtools)
devtools::install_github("https://github.com/DanielaOrschanski/OMICsdo.Stanford")
library(OMICsdo.Stanford)

### 2. Check system dependencies:
(bash)
sudo apt install cutadapt wget tar


### 3. Set up environment:
#### You will have to select the directory where all the softwares will be stored.

OMICsdo.Stanford::setup_environment()

### 4. Identify gene fusions:
RNAseqP(patients_dir = "~/Samples", 
        soft_directory = path.expand("~/OMICsdo.Stanford-Sof"), 
        cohort_name = "Example"
        )

### 5. Calculate  basic fusion stats:
fusionStats(patients_dir = "~/Samples",
            Metadata = metadata_file, 
            group = "tissue", 
            cohorte = "Example"
            )


### 6. Quality Control:
QC_seq_align(patients_dir =  "~/Elmer/Rstudio/Daniela/Muestras")

