# RiboEEF3
 Ribosome profiling analysis of eEF3-depleted Saccharomyces cerevisiae

### Prerequisites

1) Python:

  Anaconda Python v.3.9 from [Continuum](https://www.continuum.io/downloads). It comes with a bunch of libraries and have a nice  package manager `conda`. Before `conda` is able to install bioinformatic libraries/programs you have add the _bioconda_ channel.

```
    conda config --add channels bioconda
```

2) [cutadapt](https://cutadapt.readthedocs.io/en/stable/)

```
    conda install bioconda:cutadapt
```

3) `pigz` - optional if not installed cutadapt falls back to single core mode

4) [HISAT2](ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads)

```
    conda install bioconda::hisat2
```   
   Hisat2 trims ends of reads with bad quality by default. That leads to uncorrect mapping of ribosome location. From the version 2.0.5 there is an option to turn this behavior off.

5) [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

```
    conda install bioconda::bowtie2
```   

6) [samtools](https://github.com/samtools/samtools/)

```
    conda install bioconda::samtools
```

7) [pysam](https://github.com/pysam-developers/pysam)

```
    conda install bioconda::pysam
```

### Install dependencies using pre-created conda env config

  An easier way to replicate the packages and its dependencies is to just use the environment.yml to create an identical conda environment.

```
    conda env create -f environment.yml
``` 

### Running the scripts

  To run the script in the CREATE HPC, execute the following steps.

1) Clone the latest version of the repo

```
    git clone https://github.com/bpulugundla/RiboEEF3.git
```

2) Go to the directory with the repo and then create a run script with 

```
    #!/bin/bash

    #SBATCH --job-name=RNASeq 
    #SBATCH --output=rnaseq.out
    #SBATCH --time=5:00:00
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=6
    #SBATCH --mem-per-cpu=1GB

    cd RiboEEF3

    srun ./run.sh --refdir References/ --rawdir /scratch_tmp/grp/msc_appbio/group3/data/Raw/ --groupA "ERR2660269 ERR2660271" --groupB "ERR2660264 ERR2660265" --stage 2 --stop_stage 7
```

Note: The start stage for the RNA-Seq samples is 2, since we skip the adapter trimming stage. Group A points to the depleted samples and group B to the wild type.

3) Run the above created bash script with `sbatch -p msc_appbio <the_above_run_script>` from outside the repo directory.

4) The output files are stored in `/scratch_tmp/grp/msc_appbio/group3/data/`. The logs at each stage are in its respective `Reports` sub directory.

### Additional data files

  * Genome.fa  - genome sequence in FastA format
  * ncRNA.fa   - non coding RNA in FastA format
  * Genome.gtf - genome annotation in GTF (gff2) format

Other data files are derived based on those three and commands for that are described in the file  `build_index.sh`.
_Saccharomyces cerevisiae_ genome, annotation, ncRNA and indexes are locating in the folder **References/**.
