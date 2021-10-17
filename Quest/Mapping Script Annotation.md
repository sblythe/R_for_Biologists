# Mapping Script Annotation

We have provided you with some test scripts (and test data) that you can use to verify that the Quest environment is set up properly, and that you know how to run a job. These scripts can easily be modified to accommodate any mapping job, however, it is important to understand how each of the parts work.

In this example, we will go through section-by-section the _paired-end_ mapping script. For those of you interested in single-end data only, this is still useful. Regardless of what you need, I encourage you to compare both the single- and paired-end scripts as you familiarize yourself with this process.

## Script setup

The top of the script has information that sets up how the entire script will be run.

```{}
#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=00:20:00
#SBATCH --mem-per-cpu=12G
#SBATCH --job-name=test_script_PE     ## change this to job name
#SBATCH --output=outlog               ## will write files to the current dir.
#SBATCH --error=errlog
```

In this section we tell the computer that what follows is to be interpreted using `bash` by giving the command `#!/bin/bash`, of which `/bin/bash` is the path to `bash`. This is how all bash scripts are initialized.

The lines that follow beginning with `#SBATCH` are used for the 'slurm' job manager on Quest. If you were running this locally on your own computer, or independently of the job scheduler, you would eliminate all these lines.

I refer you to the Quest Genomics documentation to parse what these lines mean in general. However, I will point out what lines you should edit as you are generating reads for this class.

The line `#SBATCH --job-name=test_script_PE` should be edited so that you give any new jobs a unique name. Not the end of the world if you forget to do this, but it is tidy to do so.

The two lines that follow, which specify where to write out the output and the error are very useful in that they log the two `bash` output streams, `stdout` and `stderr`. The software we use to trim and map reads routinely output information to both of these streams (and not all `stderr` output is necessarily an error). These are good to capture, and the way they are specified here is that they will be written to a file named either `outlog` or `errlog` that will be located in the current directory at the time of script execution. This is, in other words, sloppy. In practice, you can specify an entire filepath after `#SBATCH --output=` so that you have control over where these files are written.

## Specifying the variables

```{}
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
```

Once you have tailored this code to work for your particular files, so long as you are consistent in the way you name files (more on this later), these lines here should be the only ones that you need to check and change as you process different and new datasets.

There are four things going on in this chunk. 1) We assign to the variable `basedir` the path to the base directory where our raw and processed data will go. You have to change this to the full path to that directory. 2) We make new directories within the base directory that will hold the output from trimming and mapping. 3) We load the necessary modules. 4) We assign to variables the remaining necessary objects for specifying input and output within the code that follows. As noted in the comment, the last object, `index`, points to a location where I have generated a bowtie index for the _Drosophila_ genome version dm6. If you are mapping to the mouse genome, you want to find the location of those indices and specify the correct location here. Note the syntax of this line is `path/to/index/common-suffix-of-the-6-index-files`.

Implicit in all of this is the directory structure that I suggest you use for this process.

1) Each experiment gets its own base directory.

2) The raw .fastq files for the experiment resides in the base directory within a subdirectory named `Raw_Data`.

3) The trimmed and mapped data will be in their own subdirectories, `Trimmed_Reads` and `Mapped_Reads`.

4) For future reference, the script you run to generate these data **should be stored within the base directory as well**.

As currently written, if you set up a base directory, copy the raw data to the subdirectory `Raw_Data`, and place the script in the base directory, then you should have at the end of the script:

```{}
>my_base_directory
    >Raw_Data
    >Trimmed_Reads
    >Mapped_Reads
    -outlog
    -errlog
    -script
```

This then constitutes the sum total of your input data for downstream analysis, plus the script that generated the data, as well as the entire output of the script in the form of the output and error logs. We also generate some additional files during the analysis that will come in handy when two years down the line you need to remember how you made the data and what the mapping rate was, et cetera.

## Trimming the raw data

```{bash, eval = FALSE}
source /projects/b1059/RforBiologists/TrimGaloreEnv/bin/activate 

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
```

This chunk represents a for-loop that will take some set of files in the `Raw_Data` directory and loop through them to perform adapter trimming with `TrimGalore`. I have installed `TrimGalore` locally in the class allocation and instructions on how to do this will be provided separately. There are also alternative trimmer programs that can be used instead. I like `TrimGalore` because it is simple and quick.

`TrimGalore` runs in a virtual environment. The first line of this chunk initializes the virtual environment:

`source /projects/b1059/RforBiologists/TrimGaloreEnv/bin/activate`

Before we discuss the for-loop, let me describe one of the issues with looping over paired-end data. Unlike single-end data, paired-end data has two fastq files per sample. One corresponds to "read1" and the other corresponds to "read2". We need to present both of these files together to `TrimGalore` (and also to `bowtie2`, below). The way to set up a loop to process *both* of these files is to loop over the 'read1' files, and to specify the corresponding 'read2' file within each iteration of the for-loop.

**If you are solely doing single-end data, please compare the two scripts. While you don't have to worry about regenerating the filepath to the second read, you do need to craft the for-loop call carefully as described below.**

For-loops in `bash` have the format:

```{}
for x in *some range*
    do
        execute these commands
    done
```

We therefore want to take each file in `Raw_Data` corresponding to 'read1' and assign it to a variable that we can manipulate within the for-loop. The test script is set up to deal with *the test data that I have provided*. When applying this script to *your* data, you will want to pay special attention to the following line that initializes the for-loop.

```{}
for file in ${basedir}/${rawdir}/*_R1_2M.fastq.gz
```

In plain language, what the line above says is, for each file in the raw data directory located in the base directory that ends with the text "_R1_2M.fastq.gz" (`${basedir}/${rawdir}/*_R1_2M.fastq.gz`), assign that filepath to the object `file` and iterate through those in the following loop.

When you have your own files, you want to examine all the 'read1' files and determine what generalizations you can apply to list all of the 'read1' files and not any 'read2' files. To experiment with this, you can use the `ls` function in terminal.

Another note is the `*` character, which in `bash` stands for a 'wildcard'. The `*` stands in for as few as zero additional characters of any kind. It can appear at the beginning, in the middle, or at the end of a character string to specify a filename.

Moving on: within the for-loop (i.e., within `do` and `done`), the steps are that we assign to `filename` the base name of the file (i.e., we extract the filename from the full path). We do that with this command:

```{}
filename=$(basename "$file")
```

We then report back to `stdout` some sanity checks that tell us that we can properly specify the filepath to the file of interest using solely the variables we have assigned. To do this, we use the `bash` operation `echo`, which just prints to `stdout`.

Next, we generate the filename of 'read2' from that of 'read1' by substituting any text in the character string that differs between the two. This is done with this line:

```{}
filename2=$(echo $filename | sed "s/R1/R2/g")
```

In plain language, we are assigning to `filename2` the result of a piped operation. First we `echo` `$filename`, and `echo` prints to `stdout`. However, we now invoke the `pipe` (`|`), which passes `stdout` and makes it the input (`stdin`) to the subsequent command, which in this case is `sed` (= string edit). `sed` takes the filename and exchanges the characters "R1" for "R2" using its expected syntax. Please see the manual for `sed` for more details (`man sed`).

We then do another set of `echo`s to print to `stdout` additional checks that we can indeed specify the correct path to file 2.

**One way that these scripts fail is that you mess up either the for-loop specification, or the formation of the filename for file 2. This information will be saved in the `outlog` file and you can refer to that if you run into problems.**

Next, we finally come to running `TrimGalore`:

```{}
trim_galore \
   -o ${basedir}/${trimdir} \
   --gzip \
   --fastqc \
   --cores 8 \
   --paired ${basedir}/${rawdir}/${filename} ${basedir}/${rawdir}/${filename2}
```

This is one command that I have broken up into multiple lines for readability. This is generally good practice for making readable code. To do this, I use the backslash (`\`). Here we have our first call to a `bash` function. The first command (`trim_galore`) calls the function we want to use. The rest of the command provides options and required input to the function.

In plain language: we are running TrimGalore (`trim_galore`) and whatever output it generates will be sent to the indicated output directory (`-o ${basedir}/${trimdir}`). We will be compressing the output (`--gzip`), running fastQC on the trimmed data (`--fastqc`), and doing all of this with 8 (the maximum suggested) cores/threads (`--cores 8`). The input data will be paired end, and we specify the location of the two input files. `--paired ${basedir}/${rawdir}/${filename} ${basedir}/${rawdir}/${filename2}`.

To explore the options available for running `trim_galore`, first activate the virtual environment as described above, and then enter `trim_galore` or `trim_galore -h` in the command line. It will print a detailed description of options to the screen.

The loop will loop through all files you specify, generating trimmed, validated paired-end reads. With default settings, `TrimGalore` will drop both pairs if one of them is trimmed to less than 20 bp in length. It will automatically identify any contaminating Illumina sequencing adapter and remove it. You can specify specific ones if you like (see options).

Once the looping is done, we no longer need the virtual environment, so we stop it:

```{}
deactivate
```

## Mapping reads with Bowtie2

The next for-loop finishes off our mapping process, and it accomplishes a number of tasks:

1) It maps paired-end reads using `bowtie2`, outputting `.sam`.

    * It also saves the mapping rate message `bowtie2` passes to `stderr`.

2) It converts the `.sam` output to `.bam` using `samtools view`.

3) It sorts the reads according to map position using `samtools sort`.

4) It cleans up the resulting sorted `.bam` file using `Picard CleanSam`.

5) It marks PCR/Optical duplicates (but does not remove them) using `Picard MarkDuplicates`.

6) It generates an index using `samtools index`.

7) It generates a summary of the "sam flags" using `samtools flagstat`.

8) Finally, it discards intermediate files and yields, for each input pair of trimmed fastq files:

    * a `.bam` file of mapped, sorted, duplicate marked reads (`*.mapped.md.bam`)

    * the `bowtie2` mapping log.

    * the duplication metrics file from `Picard`.

    * the flagstat summary text file

    * the `.bam` index (`*.mapped.md.bam.bai`)

Here is the code:

```{}
for file in ${basedir}/${trimdir}/*R1_2M_val_1.fq.gz
  do
    filename=$(basename "$file")
    filename2=$(echo $filename |sed "s/R1/R2/g")
    filename2=$(echo $filename2 |sed "s/val_1/val_2/g")

    echo "Next, we test Bowtie2 on the output of TrimGalore:"
    echo ${basedir}/${trimdir}/${filename}
    echo ${basedir}/${trimdir}/${filename2}

    shortname=$(echo $filename |cut -d_ -f1-2)
    echo
    echo "The output for this step will have a filename that contains"
    echo ${shortname}

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
```

We face the same challenge as with the `TrimGalore` step in terms of formulating both the set of 'read1' as well as 'read2' filenames for the loop. The call to the for-loop takes into account the standard nomenclature that `TrimGalore` uses when outputting validated trimmed paired-end reads.

```{}
for file in ${basedir}/${trimdir}/*R1_2M_val_1.fq.gz
```

Within the loop, in addition to using an identical strategy for defining the filename for the read2 data, we do something new. All of the output files will have a unique prefix that identifies it. The input data from an Illumina run often contains unnecessary information. Here, we assign to `shortname` a trimmed character string that drops unnecessary characters.

```{}
shortname=$(echo $filename |cut -d_ -f1-2)
```

Again, we are `echo`ing the filename we have formulated, and again we `pipe` the output of this to another function. Here we are piping to `cut`, which simply takes a character string, finds a delimiter --in this case, an underscore-- (`d_`), and returns fields 1 and 2, as delimited by the selected delimiter (`-f1-2`). Please see `man cut` for additional information.

As with before **this is a step that can cause errors.** With consistent file-naming, however, this line can be written such that it always drops the extra ("R1_2M.fastq.gz") text on the input filenames. We print to `stdout` the result of this operation, you can check it in the `outlog`, and if your script runs with funny output filenames, this is the operation to troubleshoot. Note also that *`bowtie2` will output only one `.sam/.bam` file for a paired-end dataset*. No more pairs of files.

Next we have a call to `bowtie2` piped to `samtools` to do the `.sam` to `.bam` conversion in one go:

```{}
bowtie2 \
    -p 12 \
    -X 2000 \
    -x $index \
    -1 ${basedir}/${trimdir}/${filename} \
    -2 ${basedir}/${trimdir}/${filename2} \
    2> ${basedir}/${mapdir}/${shortname}.bowtie.mapping.log.txt | \
    samtools view \
    -bS - > ${basedir}/${mapdir}/${shortname}.mapped.bam
```

In plain language: we are invoking `bowtie2` to do mapping. This mapping will use 12 cores `-p 12` for this task. As we will see, these are paired-end reads, and we would like to limit the output to pairs that map witin 2000 bp of one another (-X 2000). To do mapping, we need to refer to the index located at the path we specified way at the top of this script (`-x $index`). The files corresponding to the first and second reads in this paired end dataset are located at `-1 ${basedir}/${trimdir}/${filename} -2 ${basedir}/${trimdir}/${filename2}`.

That could suffice for a call to `bowtie2` (we'd have to tell it where to write the output file to, however... *not shown here*). Instead, we do something a little different, and it involves some non linearity in what is going on. In `bash`, the `>` character means "write `stdout` to this file (followed by a filename). Preceding it with a "2", as in `2>` ,says, "write `stderr` to this file, and continue as you were with respect to outputting `stdout`". So here, we write `stderr` (which is the output stream to which `bowtie2` prints its message about mapping statistics) to a file for later reference using:

```{}
2> ${basedir}/${mapdir}/${shortname}.bowtie.mapping.log.txt
```

This is then followed by a `pipe`, which transfers `stdout` (i.e., the mapped data) to `samtools view`. In plain language: all the mapped data is coming through the pipe and is presented to `samtools view` which can be specified to take a piped input using the subtle option `-` (a hyphen). The input data is `.sam` formatted, but we want to output `.bam`, which we achieve by combining the options `-bS`, which could just as easily be written `-b -S`. Using the `>` operator, we write the output (`stdout`) to the file `${basedir}/${mapdir}/${shortname}.mapped.bam`.

This gets us our initial mapping, but the entries are in the order they were imaged on the flowcell, and not in the order of where they map on the genome. We will fix that next.

**I should say at this point that this is but one way to do your call to `bowtie2`. I encourage you to really read the `bowtie2` manual and explore some of the available options, particularly those regarding the sensitivity/effort of the mapper, what to do with multimappers, as well as whether the mapping should be end-to-end or locally anchored.**

Next, we use a call to `samtools sort` that will rearrange the reads in the order with which the chromosomes were listed in the reference genome we used for making the index.

```{}
samtools sort -@ 12 -m 2G -T ${basedir}/${mapdir}/temp \
    -o ${basedir}/${mapdir}/temp1.bam \
    ${basedir}/${mapdir}/${shortname}.mapped.bam
```

In plain language, we are using `samtools sort` to sort our reads located at `${basedir}/${mapdir}/${shortname}.mapped.bam`. We will be employing 12 cores/threads (`-@ 12`) and using up to 2 gigabytes of memory (`-m 2G`). When sorting occurs, the program sets aside bins of reads in temporary files that it erases when the task is done. Those temporary files will be named `-T ${basedir}/${mapdir}/temp`. The sorted `.bam` file will be saved at `-o ${basedir}/${mapdir}/temp1.bam`. This might seem a little confusing: I have chosen to output the sorted reads to a file named `temp1.bam`. The reason for this is that we are still going to be changing things about this `.bam` file before it is done, and we don't want to keep this intermediate. Hence the name. We *could have* piped the output of the prior step (`samtools view`) into this one, but I have found sometimes that I don't want that (like, I do a step between `view` and `sort`, and so historically I haven't had the `pipe` between these two steps.)

We can't `pipe` the next steps however. We now go to the `Picard` program to run `CleanSam` followed by `MarkDuplicates`. `Picard` is a `java` program and runs with slightly different syntax. As far as I can tell there is no way to easily pipe `stdout` to a call to `java`, but I could be wrong. Here's the code:

```{}
java -jar /software/picard/2.21.4/bin/picard.jar CleanSam \
    INPUT=${basedir}/${mapdir}/temp1.bam \
    OUTPUT=${basedir}/${mapdir}/temp2.bam

    java -jar /software/picard/2.21.4/bin/picard.jar MarkDuplicates \
    INPUT=${basedir}/${mapdir}/temp2.bam \
    OUTPUT=${basedir}/${mapdir}/${shortname}.mapped.md.bam \
    METRICS_FILE=${basedir}/${mapdir}/${shortname}.dup.metrics.txt
```

In plain language, we first invoke `Picard CleanSam` (`java -jar /software/picard/2.21.4/bin/picard.jar MarkDuplicates`), taking as input (`INPUT=`) the file path to the output of the call to `samtools sort` (`.../temp1.bam`), and outputting a new temporary file, `OUTPUT=.../temp2.bam`. Next we invoke `Picard MarkDuplicates` (`java -jar /software/picard/2.21.4/bin/picard.jar MarkDuplicates`) taking as input the output of `CleanSam`, and outputting two files, 1) the `OUTPUT`, which we give a nice name (`${basedir}/${mapdir}/${shortname}.mapped.md.bam`), as well as a `METRICS_FILE`, which has a distinguishing filename itself.

This essentially **completes** the mapping process. However, we do three more things that will help in the long run.

Next we discard the temporary files as well as the original mapping. The duplicate-marked `.bam` file contains the same information as the original file, but has the added benefit of being sorted as well as having its duplicates marked in case we want to remove them. Here:

```{}
rm ${basedir}/${mapdir}/temp1.bam \
    ${basedir}/${mapdir}/temp2.bam \
    ${basedir}/${mapdir}/${shortname}.mapped.bam
```

Which should be self-explanatory how we remove (`rm`) these files.

Next, we use `samtools flagstat` to write a text file containing a breakdown of the `.sam flags` for our file. Please see [this site](https://broadinstitute.github.io/picard/explain-flags.html) for further information about flags.

```{}
samtools flagstat ${basedir}/${mapdir}/${shortname}.mapped.md.bam > \
    ${basedir}/${mapdir}/${shortname}.mapped.md.bam.FLAGSTAT.txt
```

And finally, we create a `.bam index` for our file. These indices are not obviously useful, but many downstream applications use them (and expect to find them in the same directory as the `.bam` file) when accessing the data. It is supremely annoying to have to stop an analysis to create indices. Let's do it here.

```{}
samtools index ${basedir}/${mapdir}/${shortname}.mapped.md.bam
```

This completes one round of the loop.

## Afterword

You are putting forth a lot of effort at this moment just to do the first step of your genomics analysis. By the time all is said and done, trimming and mapping will comprise just 1% of the time you spend with this data, and these steps will be a distant memory. You will probably never touch the raw or the trimmed reads ever again. However, once you go to publish these data, all of a sudden you need to have access to these files again, and not only that-- you need to provide numbers.

You *should* report in a publication the total number of reads per sample, mapping rates, and the number of reads that are discarded due to quality issues or trimming. You *have to* report these numbers when you submit the raw data to Gene Expression Omnibus/Sequence Read Archive. Nothing is more depressing than thinking you are a day away from submitting a manuscript than realizing that you first have to get a GEO accession number, and then seeing how much work it is to put together the submission to GEO.

In the example script, I have included some calls to `echo` that will report the version of all the critical programs we just used. These will be recorded in `outlog` and may not seem valuable now, but will be incredibly valuable when you are writing your incredibly detailed Materials and Methods section of your first paper/dissertation. This information is also requested by GEO in cases where you submit processed data. Please remember to give credit where credit is due. Each of these applications are cite-able, and you should be sure to do so when you write up your results.

The critical numbers are present in these output files that we have generated. Total number of reads in the original and trimmed fastQ files can be found in the `TrimGalore` trimming report, located in the `Trimmed_Reads` directory. The overall mapping rate, as well as more nuanced breakdowns of how things mapped can be found in the `Mapped_Reads` directory either in the `mapping.log.txt` or the `flagstat` output. In a worst-case scenario, you are cutting and pasting these numbers into the GEO Submission spreadsheet. In a best case scenario, you learn a ton of `bash` scripting and can write a function that will scrape this information from these files and populate a tab-delimited text document. In whatever case, having this data in a predictable location now is waaaay better than needing to figure out at the last minute how to generate these metrics when you are months or years removed from the initial mapping process.

In any case, the way that this script is set up, it should not only help move your analysis along at this early stage of the research, but should also help reduce headaches at the end of the project. This is definitely something that is missed when the mapping is done entirely in Galaxy, or if the mapping is done in an improvised manner in `bash` without keeping all the relevant files and information.
