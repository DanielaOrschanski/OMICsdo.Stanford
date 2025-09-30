# OMICsdo.Stanford

R package for RNAseq data analysis.

Using the command OMICSdo.Stanford::setup_environment() it downloads the following softwares:
1. FastQC
2. Trim Galore
3. STAR
4. Arriba
5. R packages required

Additionally, it downloads the hg38 genome reference (fasta and annotation).

The only "complicated" step is the index of the reference with STAR. In order to do this we provide 2 options:
1. Downloading it from Google Drive Folder: setting parameter "fromDrive = TRUE" in downloadHG38() function.
2. Letting STAR generates its index. No need for changes in any parameter.

The system packages (installed by sudo) required are:
sudo apt install wget unzip cutadapt tar bzip2

