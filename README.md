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

### Additional data files

  * Genome.fa  - genome sequence in FastA format
  * ncRNA.fa   - non coding RNA in FastA format
  * Genome.gtf - genome annotation in GTF (gff2) format

Other data files are derived based on those three and commands for that are described in the file  `build_index.sh`.
_Saccharomyces cerevisiae_ genome, annotation, ncRNA and indexes are locating in the folder **References/**.
