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

