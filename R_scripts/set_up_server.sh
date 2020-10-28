#'#################################################################################
#'#################################################################################
#' Set up server for MOI project
#'#################################################################################
#'#################################################################################

## Install HaplotypR (R) ####
BiocManager::install("ShortRead")
devtools::install_github("lerch-a/Rswarm")
devtools::install_github("lerch-a/Rvsearch")
devtools::install_github("lerch-a/HaplotypR")

## Install Miniconda3
cd software/
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 766 Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh 
# Follow instructions
### Install in /home/isglobal.lan/cruiz/miniconda3
## Add links to bin folder
ln -s ~/miniconda3/bin/conda ~/bin/conda
ln -s ~/miniconda3/bin/python3.7 ~/bin/python3
conda install -c conda-forge biopython ## Install biopython (for other dependencies)

## Install seekDeep (home)
cd software/
git clone https://github.com/bailey-lab/SeekDeep
cd SeekDeep
./configure.py 
./setup.py --compfile compfile.mk --outMakefile makefile-common.mk
make

## Install bowtie2
conda install -c bioconda bowtie2 
## Add link to bin folder
ln -s ~/miniconda3/bin/bowtie2 ~/bin/bowtie2
ln -s ~/miniconda3/bin/bowtie2-build ~/bin/bowtie2-build

## Download phiX174 reference in FASTA (for bowtie2)
mkdir ~/PublicData/REFERENCES/phiX174/
cd ~/PublicData/REFERENCES/phiX174/
nano phiX174.fa
### Copy content from https://www.ncbi.nlm.nih.gov/nuccore/J02482.1?report=fasta

### Create bowtie2 index
bowtie2-build phiX174.fa phiX174


# Download plasmodium references
mkdir ~/PublicData/REFERENCES/PlasmodiumFalciparum/
cd ~/PublicData/REFERENCES/PlasmodiumFalciparum/
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/DD2/Assembly/V1_morphed/PfDD2_morphed_v1.union.embl.gz
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/HB3/Assembly/V1_morphed/PfHB3_morphed_CV.union.embl.gz
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/7G8/Assembly/V1_morphed/Pf7G8_morphed_v1.union.embl.gz
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/3D7/3D7.latest_version/version3.1/2019/January_2019/*.*

## Convert to FASTA to allow indexing with bowtie
### Decompres
for i in `ls *.gz`
do
 gzip -d $i
done

### Use the following script (https://github.com/etal/biofrills/blob/master/scripts/seq-convert)
nano ~/software/GenomeSequenceConverter.py
for i in `ls *.embl`
do
   name=${i%.embl}
   python ~/software/GenomeSequenceConverter.py embl fasta < $name.embl > $name.fasta
done

## Create bowtie2 index
### 3D7
bowtie2-build `ls Pf3D7*.fasta | paste -sd "," -` Pf3D7

### 7G8
bowtie2-build Pf7G8_morphed_v1.union.fasta Pf7G8

### DD2
bowtie2-build PfDD2_morphed_v1.union.fasta PfDD2

### HB3
bowtie2-build PfHB3_morphed_CV.union.fasta PfHB3

## Download raw data
mkdir ~/data/malaria
mkdir ~/data/malaria/MOI
mkdir ~/data/malaria/MOI/data
cd ~/data/malaria/MOI/data

## Add user and password from LastPass
wget -r -nd --no-parent http://seq.crg.es/download/external/Miro/Alfredo_Mayor/Himanshu_Gupta/2019-04-22_MiSeq/ --http-user=## --http-password=##
wget -r -nd --no-parent http://seq.crg.es/download/external/Miro/Alfredo_Mayor/Himanshu_Gupta/2019-06-29_MiSeq/ --http-user=## --http-password=##
wget -r -nd --no-parent http://seq.crg.es/download/external/Miro/Alfredo_Mayor/Himanshu_Gupta/2019-07-13_MiSeq/ --http-user=## --http-password=##

## Remove html files 
rm *html*

## Create file with list of files per run
wget --spider -r --no-parent http://seq.crg.es/download/external/Miro/Alfredo_Mayor/Himanshu_Gupta/2019-04-22_MiSeq/ --http-user=## --http-password=## &> run1_raw.txt
wget --spider -r --no-parent http://seq.crg.es/download/external/Miro/Alfredo_Mayor/Himanshu_Gupta/2019-06-29_MiSeq/ --http-user=## --http-password=## &> run2_raw.txt
wget --spider -r --no-parent http://seq.crg.es/download/external/Miro/Alfredo_Mayor/Himanshu_Gupta/2019-07-13_MiSeq/ --http-user=## --http-password=## &> run3_raw.txt

grep 'http.*gz$' run1_raw.txt | sed 's/.*http.*MiSeq\///' > run1.txt
grep 'http.*gz$' run2_raw.txt | sed 's/.*http.*MiSeq\///' > run2.txt
grep 'http.*gz$' run3_raw.txt | sed 's/.*http.*MiSeq\///' > run3.txt


# Install fastQC
cd ~/software
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
unzip fastqc_v0.11.8.zip
chmod 755 ./FastQC/fastqc
## Add link to bin folder
ln -s ~/software/FastQC/fastqc ~/bin/fastqc

# Install bbmap
cd software
wget -O BBMap_38.59.tar.gz https://downloads.sourceforge.net/project/bbmap/BBMap_38.59.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbbmap%2Ffiles%2Flatest%2Fdownload&ts=1563805170
tar -xvzf BBMap_38.59.tar.gz
 
# Create scratch folder
mkdir /scratch/cruiz/32_19_MOI
ln -s ~/data/malaria/MOI/data /scratch/cruiz/32_19_MOI/data 
cd  /scratch/cruiz/32_19_MOI
mkdir results
mkdir results/QC
mkdir results/HaplotypR_config
mkdir results/HaplotypR_files
mkdir results/Dillution
mkdir results/SEVMAL
mkdir results/MDA

# Add win winScp MOI_markerFile.txt to results/HaplotypR_files. This file contains 
# the primer sequence and the gene references
nano results/HaplotypR_files/MOI_markerFile.txt 
## This step ensure that the file complies with UNIX format


