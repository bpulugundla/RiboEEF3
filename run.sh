#!/usr/bin/env bash

eval "$(conda shell.bash hook)"
conda activate appbio_project

log() {
    local fname=${BASH_SOURCE[1]##*/}
    echo -e "$(date '+%Y-%m-%dT%H:%M:%S') (${fname}:${BASH_LINENO[0]}:${FUNCNAME[1]}) $*"
}

# General configuration
stage=1                 # Processes starts from the specified stage.
stop_stage=0            # Processes is stopped at the specified stage.
skip_align=false        # Skip alignment stages.
skip_eval=false         # Skip evaluation stages.
ncpu=6                  # The number of cpus.
python=python3          # Specify python to execute espnet commands.
clean=gzip              # Cleaning processed files.
quality=0.80            # Quality filter threshold for retaining reads. 
mapping=5               # Mapping reads
read_len_min=25         # Minimum length of reads included
read_len_max=35         # Maximum length of reads included
normalised=rpm          # Read data type
mapped_twice=false      # Reads mapped twice
metag_threshold=30      # Raw counts around start or stop if less ignore a gene/region
metag_span=60           # Number of nucleotides before and after start/stop in metagene profiles
metag_spancorr=90       # number of nucleotides before and after start/stop in metagenomic profile from corrected data
gene_raw_mean_thr=0.3   # Percentage of raw counts per 100 nt
codon_raw_mean_thr=1.6  # Raw counts per codon. 1.6 indicates 6 counts
codons_before_psite=1   #

# Data related
refdir=References       # Directory of reference genome.
rawdir=Raw              # Directory of raw FASTQ files.
groupA=
groupB=

help_message=$(cat << EOF
Usage: $0 --reference "<reference_dir_name>" --raw "<raw_dir_name> --groupA "MetS2 MetS4" --groupB "WTS1 WTS3"

Options:
    # General configuration
    --stage             # Processes starts from the specified stage (default="${stage}").
    --stop_stage        # Processes is stopped at the specified stage (default="${stop_stage}").
    --skip_align        # Skip alignment stages (default="${skip_align}").
    --skip_eval         # Skip evaluation stages (default="${skip_eval}").
    --ncpu              # The number of cpus (default="${ncpu}").
    --python            # Specify python version (default="${python}").
    --clean             # bgzip - Gzip files after they are no longer necessary using multiple threads - 6 is hard coded
                        # gzip  - Gzip files after they are no longer necessary
                        # rm - Remove files after they are no longer necessary (default="${clean}").  
    --quality           # Quality filter threshold for retaining reads based on PHRED scores. Higher is more selective (default="${quality}").
    --mapping           # 5 maps reads according to their 5' end. 3 maps reads according to their 3' end (default="${mapping}").
    --read_len_min      # Minimum length of reads included (default="${read_len_min}").
    --read_len_max      # Maximum length of reads included (default="${read_len_max}").
    --normalised        # rpm - use normalized data (reads per million raw), raw - raw counts (default="${normalised}").
    --mapped_twice      # Reads mapped once or twice (default="${mapped_twice}").
    --metag_threshold   # Raw counts around start or stop if less ignore a gene/region (default="${metag_threshold}").
                        # If Normalized = "rpm"  MetagThreshold is readjusted by dividing normalization_factor inside program FILTER 3.
    --metag_span        # Number of nucleotides before and after start/stop in metagene profiles (default="${metag_span}").
    --metag_spancorr    # Number of nucleotides before and after start/stop in metagenomic profile from corrected data 
                        # (default="${metag_span}").
    --gene_raw_mean_thr # Percentage of raw counts per 100 nt. - low threshold (default="${gene_raw_mean_thr}").
                        # FILTER 1  - is used in functions codonTablesA & codonTablesB converted to GeneRpmMeanThr
                        # GeneRpmMeanThr = GeneRawMeanThr / norm_factor
    --codon_raw_mean_thr# Raw counts per codon. 1.666 = 6 raw counts per codon (default="${codon_raw_mean_thr}").
                        # FILTER 2  - is used in functions codonTablesA & codonTablesB converted to CodonRpmMeanThr
                        # CodonRpmMeanThr = CodonRawMeanThr / norm_factor
    --codons_before_psite # affects sequence part of codon tables and Master table. 1 - A-site included (default="${codons_before_psite}").
                        # 2 - A-site + one codon before

    # Data related
    --refdir            # Directory of reference genome (default="${refdir}").
    --rawdir            # Directory of raw FASTQ files (default="${rawdir}").
    --groupA            # Group refers to samples of different replicas but same conditions
                        # Replica is defined as order of sample names under Group
                        # Space separated samples incubated under specific conditions
    --groupB

EOF
)

log "$0 $*"
# Save command line args for logging (they will be lost after utils/parse_options.sh)
run_args=$(utils/print_args.sh $0 "$@")
. utils/parse_options.sh

echo "${run_args}"

basedir=$(dirname "${rawdir}")
samples="${groupA} ${groupB}" 

[ -z "${groupA}" ] && echo "Error: --groupA samples are required." && exit 1;
[ -z "${groupB}" ] && echo "Error: --groupB samples are required." && exit 1;

for name in ${groupA}; do
    [ ! -f ${rawdir}/$name.fastq.gz ] && echo "$0: Sample ${name} not found at ${rawdir}/$name" && exit 1;
done

for name in ${groupB}; do
    [ ! -f ${rawdir}/$name.fastq.gz ] && echo "$0: Sample ${name} not found at ${rawdir}/$name" && exit 1;
done

if [ ${stage} -le 1 ] && [ ${stop_stage} -ge 1 ]; then
    log "Step 1: Trimming adapter sequnces from the ${samples} FASTQ files in ${rawdir} using CutAdapt. CutAdapt removes adapter sequences from high-throughput sequencing reads."
    mkdir -p ${basedir}/1-Trimmed ${basedir}/1-Trimmed/Reports
    for sample in ${samples}; do
        input="${rawdir}/${sample}.fastq.gz"
        output="${basedir}/1-Trimmed/${sample}_trimmed.fastq.gz"
        cutadapt --discard-untrimmed -j ${ncpu} -a "CTGTAGGCACCATCAAT" \
            -o ${output} ${input} > ${basedir}/1-Trimmed/Reports/${sample}.txt || exit 1
    done
fi

if [ ${stage} -le 2 ] && [ ${stop_stage} -ge 2 ]; then
    log "Step 2: Quality filter the trimmed FASTQ files based on the PHRED scores."
    ${python} local/qc.py --rawdir ${rawdir} --names "${samples}" --threshold ${quality} \
        --read_len_min ${read_len_min} --read_len_max ${read_len_max}
fi

if [ ${stage} -le 3 ] && [ ${stop_stage} -ge 3 ]; then
    log "Step 3: Remove non-coding RNAs sequences from the filtered FASTQ files using Bowtie2."
    mkdir -p ${basedir}/3-Subtracted ${basedir}/3-Subtracted/Reports ${basedir}/3-Subtracted/SAM
    for sample in ${samples}; do 
        ncRNA="${refdir}/Indexes/ncRNA"
        input="${basedir}/2-Filtered/${sample}_filtered.fastq"
        output="${basedir}/3-Subtracted/SAM/${sample}.sam"
        unmapped="${basedir}/3-Subtracted/${sample}.fastq"
        bowtie2 --no-unal -p ${ncpu} --un ${unmapped} -x ${ncRNA} \
            -U ${input} -S ${output} > ${basedir}/3-Subtracted/Reports/${sample}.txt || exit 1
    done
fi

if [ ${stage} -le 4 ] && [ ${stop_stage} -ge 4 ]; then
    log "Step 4: Align sequences to the genome from the ncRNA subtracted FASTQ files"
    mkdir -p ${basedir}/4-Aligned ${basedir}/4-Aligned/Reports
    for sample in ${samples}; do
        genome="${refdir}/Indexes/Genome"
        input="${basedir}/3-Subtracted/${sample}.fastq"
        output="${basedir}/4-Aligned/${sample}.sam"
        hisat2 --no-unal -p ${ncpu} -k 2 --no-softclip --dta -x ${genome} \
            -U ${input} -S ${output} > ${basedir}/4-Aligned/Reports/${sample}.txt || exit 1
    done
fi
