#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=00:20:00
#SBATCH --mem-per-cpu=12G
#SBATCH --job-name=test_script_PE                 ## change this to job name
#SBATCH --output=outlog            ## will write files to the current dir.               
#SBATCH --error=errlog

basedir=</full/path/to/the/project/directory>     ## change to full path

## make destinations                 
mkdir ${basedir}/Trimmed_Reads             
mkdir ${basedir}/Mapped_Reads               

## load modules be verbose about versions since the defaults can change!
module load python/3.8.4
module load fastqc/0.11.5
module load bowtie2/2.4.1
module load samtools/1.10.1
module load picard/2.21.4

## assign the variables for the filepaths for this mapping
rawdir=Raw_Data
trimdir=Trimmed_Reads
mapdir=Mapped_Reads
index=/projects/b1059/RforBiologists/Bowtie_Indices/dm6/dm6 ## make sure you have correct index here

################################################################################
################################################################################
## ---PAIRED END DEMO--- Can be modified for routine use                      ##
##                                                                            ##
## The idea is that you would have a script where you only need to change     ##
## the variables above this line, and the stuff below always runs perfectly.  ##
## This takes some trial and error. Also, the code below is *not* generalized ##
## to handle whatever you throw at it, particularly:                          ##
## 1) the file identification at the head of the two for-loops below.         ##
## 2) the verbosity of the code (most of the 'echo' lines below). You can     ##
##    comment those lines out and retain the functional annotations.          ##
## 3) the call to `cut` in each for-loop that generates the short file name.  ##
## 4) (in your home lab allocations) you will need to change the path to the  ##
##    TrimGalore virtual environment. This should work fine for class though. ##
## ... so this is not plug and play, but if you modify those lines, it should ##
## work for you all the time.                                                 ##
################################################################################
################################################################################

## activate the virtual environment for TrimGalore/cutadapt
source /projects/b1059/RforBiologists/TrimGaloreEnv/bin/activate

## loop over files in the raw data directory that end in *.fastq.gz to trim
## adapter sequences, and perform fastQC. The big difference between single end
## and paired-end mappings is that you need to loop through the Read1 files and
## construct the filename of Read2. Both files are presented to TrimGalore. The trick
## is to reliably designate the filenames without errors. 

for file in ${basedir}/${rawdir}/*_R1_2M.fastq.gz
  do
   filename=$(basename "$file")
   echo "path to the first read"
   echo $basedir/$rawdir/$filename
   echo "base name for the first read"
   echo $filename
   filename2=$(echo $filename | sed "s/R1/R2/g") ## swaps "R1" for "R2" in filename
   echo "path to the second read"
   echo $basedir/$rawdir/$filename2
   echo "base name for the second read"
   echo $filename2
   echo
   echo

   trim_galore \
   -o ${basedir}/${trimdir} \
   --gzip \
   --fastqc \
   --cores 8 \
   --paired ${basedir}/${rawdir}/${filename} ${basedir}/${rawdir}/${filename2}

  done

## deactivate the virtual environment...
deactivate

## loop over files in the trimmed data directory that end in *R1_2M_val_1.fq.gz
## to map and sort and mark duplicates. Note that TrimGalore has a different filename
## format for paired-end reads than it does for single end. This is reflected here.
## Bowtie2 has different parameters as well for paired-end reads. 

for file in ${basedir}/${trimdir}/*R1_2M_val_1.fq.gz
  do
    filename=$(basename "$file")
    filename2=$(echo $filename |sed "s/R1/R2/g")
    filename2=$(echo $filename2 |sed "s/val_1/val_2/g")
    echo
    echo "Next, we test Bowtie2 on the output of TrimGalore:"
    echo ${basedir}/${trimdir}/${filename}
    echo ${basedir}/${trimdir}/${filename2}

    shortname=$(echo $filename |cut -d_ -f1-2)
    echo
    echo "The output for this step will have a filename that contains"
    echo ${shortname}
    echo

    bowtie2 \
    -p 12 \
    -X 2000 \
    -x $index \
    -1 ${basedir}/${trimdir}/${filename} \
    -2 ${basedir}/${trimdir}/${filename2} \
    2> ${basedir}/${mapdir}/${shortname}.bowtie.mapping.log.txt | \
    samtools view \
    -bS - > ${basedir}/${mapdir}/${shortname}.mapped.bam

    samtools sort -@ 12 -m 2G -T ${basedir}/${mapdir}/temp \
    -o ${basedir}/${mapdir}/temp1.bam \
    ${basedir}/${mapdir}/${shortname}.mapped.bam

    java -jar /software/picard/2.21.4/bin/picard.jar CleanSam \
    INPUT=${basedir}/${mapdir}/temp1.bam \
    OUTPUT=${basedir}/${mapdir}/temp2.bam

    java -jar /software/picard/2.21.4/bin/picard.jar MarkDuplicates \
    INPUT=${basedir}/${mapdir}/temp2.bam \
    OUTPUT=${basedir}/${mapdir}/${shortname}.mapped.md.bam \
    METRICS_FILE=${basedir}/${mapdir}/${shortname}.dup.metrics.txt

    rm ${basedir}/${mapdir}/temp1.bam \
    ${basedir}/${mapdir}/temp2.bam \
    ${basedir}/${mapdir}/${shortname}.mapped.bam

    samtools flagstat ${basedir}/${mapdir}/${shortname}.mapped.md.bam > \
    ${basedir}/${mapdir}/${shortname}.mapped.md.bam.FLAGSTAT.txt

    samtools index ${basedir}/${mapdir}/${shortname}.mapped.md.bam
done

## The lines below record the software versions you used for this analysis.

source /projects/b1059/RforBiologists/TrimGaloreEnv/bin/activate
echo "TrimGalore:"
echo $(trim_galore --version)
echo "Cutadapt:"
echo $(cutadapt --version)
deactivate

echo "Bowtie2:"
echo $(bowtie2 --version)
echo "Samtools:"
echo $(samtools --version)
echo "Picard:"
echo $(java -jar /software/picard/2.21.4/bin/picard.jar MarkDuplicates --version)
