#'#################################################################################
#'#################################################################################
#' Run QC sequencing of the data 
#'#################################################################################
#'#################################################################################

QC=results/QC
refs=~/PublicData/REFERENCES/

# phiX reads
## Map phiX reads to phiX genome using bowtie2
bowtie2 -x $refs/phiX174/phiX174 -1 data/Undetermined_S0_L001_R1_001.fastq.gz -2 data/Undetermined_S0_L001_R2_001.fastq.gz -p 30 -S $QC/phiX.sam
# 23375267 reads; of these:
#   23375267 (100.00%) were paired; of these:
#     9331641 (39.92%) aligned concordantly 0 times
#     14043626 (60.08%) aligned concordantly exactly 1 time
#     0 (0.00%) aligned concordantly >1 times
#     ----
#     9331641 pairs aligned concordantly 0 times; of these:
#       540956 (5.80%) aligned discordantly 1 time
#     ----
#     8790685 pairs aligned 0 times concordantly or discordantly; of these:
#       17581370 mates make up the pairs; of these:
#         16958302 (96.46%) aligned 0 times
#         623068 (3.54%) aligned exactly 1 time
#         0 (0.00%) aligned >1 times
# 63.73% overall alignment rate

## Get values with reformat
~/software/bbmap/reformat.sh in=$QC/phiX.sam qahist=$QC/phiX_qah.txt mhist=$QC/phiX_mh.txt

## Map phiX reads to P. Falciparum genome
bowtie2 -x $refs/PlasmodiumFalciparum/Pf3D7 -1 data/Undetermined_S0_L001_R1_001.fastq.gz -2 data/Undetermined_S0_L001_R2_001.fastq.gz -p 30 -S $QC/phiX_Pf.sam
# 23375267 reads; of these:
#   23375267 (100.00%) were paired; of these:
#     17733920 (75.87%) aligned concordantly 0 times
#     5641347 (24.13%) aligned concordantly exactly 1 time
#     0 (0.00%) aligned concordantly >1 times
#     ----
#     17733920 pairs aligned concordantly 0 times; of these:
#       273 (0.00%) aligned discordantly 1 time
#     ----
#     17733647 pairs aligned 0 times concordantly or discordantly; of these:
#       35467294 mates make up the pairs; of these:
#         35262662 (99.42%) aligned 0 times
#         204632 (0.58%) aligned exactly 1 time
#         0 (0.00%) aligned >1 times
# 24.57% overall alignment rate

## Get list of unmapped reads
samtools view -f 4 $QC/phiX_Pf.sam |  cut -f 1 | sort > $QC/phiX_Pf_unmapped.txt
samtools view -f 4 $QC/phiX.sam |  cut -f 1 | sort > $QC/phiX_unmapped.txt

wc $QC/*unmapped* ## Unmapped reads per file

comm -12 --check-order $QC/phiX_Pf_unmapped.txt $QC/phiX_unmapped.txt | wc

## Run FastQC on phiX reads
fastqc data/Undetermined* -o $QC
for i in `ls $QC/Undetermined*.zip`
do
  unzip $i -d $QC
  name=${i%.zip}
  ## Extract Quality distribution data
  awk '/Per base sequence quality/{flag=1;next}/END_MODULE/{flag=0}flag' $name/fastqc_data.txt > $name/fastqc_QCdata.txt
done


# HB3 dataset

## Get genes coordiantes in HB3 genome
### Create a fasta with gene sequences (Himunsha file)
nano $QC/HB3_genes.fa

### Map genes to reference
bowtie2 -x $refs/PlasmodiumFalciparum/PfHB3 -U $QC/HB3_genes.fa -f -S $QC/HB3_genes.sam

## Coords
### CSP: Transfer.PfHB3_01_morphed_C1.final:1808872-1809191
### AMA: Transfer.PfHB3_01_morphed_C1.final:14055486-14055843

## Map reads from mixture 1 (only HB3) to HB3 genome
### Mcsp
#### Set 1
bowtie2 -x $refs/PlasmodiumFalciparum/PfHB3 -1 data/33224_CMd1-1_F14_R1_McspSet1_S259_L001_R1_001.fastq.gz -2 data/33224_CMd1-1_F14_R1_McspSet1_S259_L001_R2_001.fastq.gz -p 30 -S $QC/HB3_Mcsp_Set1_lane1.sam
bowtie2 -x $refs/PlasmodiumFalciparum/PfHB3 -1 data/33224_CMd1-1_F14_R1_McspSet1_S438_L001_R1_001.fastq.gz -2 data/33224_CMd1-1_F14_R1_McspSet1_S438_L001_R2_001.fastq.gz -p 30 -S $QC/HB3_Mcsp_Set1_lane2.sam
bowtie2 -x $refs/PlasmodiumFalciparum/PfHB3 -1 data/33224_CMd1-1_F14_R1_McspSet1_S451_L001_R1_001.fastq.gz -2 data/33224_CMd1-1_F14_R1_McspSet1_S451_L001_R2_001.fastq.gz -p 30 -S $QC/HB3_Mcsp_Set1_lane3.sam

#### Set 2
bowtie2 -x $refs/PlasmodiumFalciparum/PfHB3 -1 data/33224_CMd1-1_F15_R1_McspSet2_S271_L001_R1_001.fastq.gz -2 data/33224_CMd1-1_F15_R1_McspSet2_S271_L001_R2_001.fastq.gz -p 30 -S $QC/HB3_Mcsp_Set2_lane1.sam
bowtie2 -x $refs/PlasmodiumFalciparum/PfHB3 -1 data/33224_CMd1-1_F15_R1_McspSet2_S450_L001_R1_001.fastq.gz -2 data/33224_CMd1-1_F15_R1_McspSet2_S450_L001_R2_001.fastq.gz -p 30 -S $QC/HB3_Mcsp_Set2_lane2.sam
bowtie2 -x $refs/PlasmodiumFalciparum/PfHB3 -1 data/33224_CMd1-1_F15_R1_McspSet2_S463_L001_R1_001.fastq.gz -2 data/33224_CMd1-1_F15_R1_McspSet2_S463_L001_R2_001.fastq.gz -p 30 -S $QC/HB3_Mcsp_Set2_lane3.sam

### Mama
#### Set 1
bowtie2 -x $refs/PlasmodiumFalciparum/PfHB3 -1 data/33225_AMd1-1_F16_R1_MamaSet1_S283_L001_R1_001.fastq.gz -2 data/33225_AMd1-1_F16_R1_MamaSet1_S283_L001_R2_001.fastq.gz -p 30 -S $QC/HB3_Mama_Set1_lane1.sam
bowtie2 -x $refs/PlasmodiumFalciparum/PfHB3 -1 data/33225_AMd1-1_F16_R1_MamaSet1_S462_L001_R1_001.fastq.gz -2 data/33225_AMd1-1_F16_R1_MamaSet1_S462_L001_R2_001.fastq.gz -p 30 -S $QC/HB3_Mama_Set1_lane2.sam
bowtie2 -x $refs/PlasmodiumFalciparum/PfHB3 -1 data/33225_AMd1-1_F16_R1_MamaSet1_S475_L001_R1_001.fastq.gz -2 data/33225_AMd1-1_F16_R1_MamaSet1_S475_L001_R2_001.fastq.gz -p 30 -S $QC/HB3_Mama_Set1_lane3.sam

#### Set 2
bowtie2 -x $refs/PlasmodiumFalciparum/PfHB3 -1 data/33225_AMd1-1_F17_R1_MamaSet2_S295_L001_R1_001.fastq.gz -2 data/33225_AMd1-1_F17_R1_MamaSet2_S295_L001_R2_001.fastq.gz -p 30 -S $QC/HB3_Mama_Set2_lane1.sam
bowtie2 -x $refs/PlasmodiumFalciparum/PfHB3 -1 data/33225_AMd1-1_F17_R1_MamaSet2_S474_L001_R1_001.fastq.gz -2 data/33225_AMd1-1_F17_R1_MamaSet2_S474_L001_R2_001.fastq.gz -p 30 -S $QC/HB3_Mama_Set2_lane2.sam
bowtie2 -x $refs/PlasmodiumFalciparum/PfHB3 -1 data/33225_AMd1-1_F17_R1_MamaSet2_S487_L001_R1_001.fastq.gz -2 data/33225_AMd1-1_F17_R1_MamaSet2_S487_L001_R2_001.fastq.gz -p 30 -S $QC/HB3_Mama_Set2_lane3.sam

## Get Mapping statistics
genes=( Mcsp Mama )
for gen in "${genes[@]}"
do
  echo $gen
  echo -e "Dataset \t All Reads \t Unmapped"
  for i in {1..2}
  do 
    for j in {1..3}
    do
    all=`samtools view -c $QC/HB3_${gen}_Set${i}_lane${j}.sam`
    unmapped=`samtools view -f4 -c $QC/HB3_${gen}_Set${i}_lane${j}.sam`
    echo -e "Set$i - Run$j\t$all\t$unmapped"
    done
  done
done

# Mcsp
# Dataset          All Reads       Unmapped
# Set1 - Run1     10494   20
# Set1 - Run2     7900    22
# Set1 - Run3     23188   68
# Set2 - Run1     6752    23
# Set2 - Run2     5356    17
# Set2 - Run3     15792   40
# Mama
# Dataset          All Reads       Unmapped
# Set1 - Run1     6634    169
# Set1 - Run2     4304    104
# Set1 - Run3     12286   338
# Set2 - Run1     10814   198
# Set2 - Run2     7910    154
# Set2 - Run3     24808   448

## Convert to bam and index
genes=( Mcsp Mama )
for gen in "${genes[@]}"
do
  for i in {1..2}
  do 
    for j in {1..3}
    do
      samtools view -S -b $QC/HB3_${gen}_Set${i}_lane${j}.sam > $QC/HB3_${gen}_Set${i}_lane${j}.bam
      samtools sort $QC/HB3_${gen}_Set${i}_lane${j}.bam -o $QC/HB3_${gen}_Set${i}_lane${j}.sorted.bam
      samtools index $QC/HB3_${gen}_Set${i}_lane${j}.sorted.bam
    done
  done
done

## Get values with reformat
genes=( Mcsp Mama )
for gen in "${genes[@]}"
do
  for i in {1..2}
  do 
    for j in {1..3}
    do
      ~/software/bbmap/reformat.sh in=$QC/HB3_${gen}_Set${i}_lane${j}.sam qahist=$QC/HB3_${gen}_Set${i}_lane${j}_qah.txt mhist=$QC/HB3_${gen}_Set${i}_lane${j}_mh.txt
    done
  done
done

## Run FastQC 
fastqc data/*CMd1-1_* -o $QC
fastqc data/*AMd1-1_* -o $QC

for i in `ls $QC/3*.zip`
do
  unzip $i -d $QC
  name=${i%.zip}
  ## Extract Quality distribution data
  awk '/Per base sequence quality/{flag=1;next}/END_MODULE/{flag=0}flag' $name/fastqc_data.txt > $name/fastqc_QCdata.txt
done