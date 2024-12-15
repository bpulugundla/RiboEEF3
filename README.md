# RiboEEF3
 The [study](https://www.nature.com/articles/s41598-019-39403-y), Ribosome profiling analysis of eEF3-depleted Saccharomyces cerevisiae, investigates the effects of depleting eukaryotic elongation factor 3 (eEF3) in Saccharomyces cerevisiae using ribosome profiling. Findings reveal that eEF3 depletion reduces translation elongation efficiency and impacts ribosomal stalling, particularly at P-site proline residues. However, it does not significantly impair ribosome recycling. This highlights eEF3's crucial, general role in yeast protein synthesis elongation.

## How to Run the Code

### 1. Clone the Repository
First, clone the repository to your local machine using Git:
```bash
git clone https://github.com/bpulugundla/RiboEEF3.git
cd RiboEEF3
```

### 2. Install Dependencies

### Setting Up the Conda Environment

Setting up the environment with the required packages/tools can be done in two ways. Both methods require Conda to be installed.

If Conda is not already installed, you can follow these steps:

1. Visit the [Miniconda Installation Guide](https://docs.conda.io/en/latest/miniconda.html) to download the Miniconda installer for your platform.
2. Download the appropriate installer for your operating system and follow the instructions.

For Linux/macOS:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

#### Option 1: Install from Conda channels 

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

4) [hisat2](ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads)

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

8) [cufflinks](https://cole-trapnell-lab.github.io/cufflinks/manual/)

```
    conda install bioconda::cufflinks
```

#### Option 2: Install dependencies using pre-created env config

  An easier way to replicate the packages and its dependencies is to just use the environment.yml to create an identical conda environment.

```
    conda env create -f environment.yml
``` 

## 3. Run the main script

```
./run.sh --refdir References --rawdir Raw --groupA "ERR2660266 ERR2660267" --groupB "ERR2660262 ERR2660263" --stage 1 --stop_stage 7
```
### Additional data files

  * Genome.fa  - genome sequence in FastA format
  * ncRNA.fa   - non coding RNA in FastA format
  * Genome.gtf - genome annotation in GTF (gff2) format

Other data files are derived based on those three and commands for that are described in the file  `build_index.sh`.
_Saccharomyces cerevisiae_ genome, annotation, ncRNA and indexes are locating in the folder **References/**.

### Get Help
For a complete list of arguments and their descriptions, run the following command:
```bash
./run.sh --help
```


