#!/usr/bin/env bash
"""
Author: Bhargav Pulugundla
Last Updated: 13-Dec-2024

This project focuses on replicating and validating the findings presented by Kasari et al. in their study, 
available at: https://www.nature.com/articles/s41598-019-39403-y.
"""

# Set bash to 'debug' mode, it will exit on :
# -e 'error', -u 'undefined variable', -o ... 'error in pipeline', -x 'print commands',
set -euo pipefail

# Activate the conda environment with the installed dependecies. Make sure to use the correct conda env.
eval "$(conda shell.bash hook)"
conda activate eef3_project

log() {
    local fname=${BASH_SOURCE[1]##*/}
    echo -e "$(date '+%Y-%m-%dT%H:%M:%S') (${fname}:${BASH_LINENO[0]}:${FUNCNAME[1]}) $*"
}

# General configuration
stage=1                 # Processes starts from the specified stage.
stop_stage=10            # Processes is stopped at the specified stage.
skip_align=false        # Skip alignment stages.
skip_eval=false         # Skip evaluation stages.
ncpu=6                  # The number of cpus.
python=python3          # Specify python to execute espnet commands.
clean=1                 # Cleaning processed files.
seq_type="ribo"         # Sequencing technique used
quality=0.80            # Quality filter threshold for retaining reads. 
mapping="5"             # Mapping read end
read_len_min=25         # Minimum length of reads included
read_len_max=35         # Maximum length of reads included
normalised="rpm"        # Read data type
mapped_twice=0          # Reads mapped twice, 0 for false and 1 for true
metagene_threshold=30   # Raw counts around start or stop if less ignore a gene/region
metagene_span=60        # Number of nucleotides before and after start/stop in metagene profiles
metagene_spancorr=90    # number of nucleotides before and after start/stop in metagenomic profile from corrected data
gene_raw_mean_thr=0.3   # Percentage of raw counts per 100 nt
codon_raw_mean_thr=1.6  # Raw counts per codon. 1.6 indicates 6 counts
codons_before_psite=1   #

# Data related
refdir=References       # Directory of reference genome.
rawdir=Raw              # Directory of raw FASTQ files.
groupA=""               # Depleted sample names. Space separated
groupB=""               # Normal sample names. Space separated.

help_message=$(cat << EOF
Usage: $0 --reference "<reference_dir_name>" --raw "<raw_dir_name> --groupA "MetS2 MetS4" --groupB "WTS1 WTS3"

Options:
    # General configuration
    --stage             # Processes starts from the specified stage (default="${stage}").
    --stop_stage        # Processes is stopped at the specified stage (default="${stop_stage}").
    --skip_align        # Skip alignment stages (default="${skip_align}").
    --skip_eval         # Skip evaluation stages (default="${skip_eval}").
    --ncpu              # The number of cpus for parallel processing (default="${ncpu}").
    --python            # Specify python version (default="${python}").
    --clean             # Remove files after they are no longer necessary (default="${clean}").
    --seq_type          # Type of sequencing used, RiboSeq(ribo) or RNA-Seq(rna) (default="${seq_type}").
    --quality           # Quality filter threshold for retaining reads based on PHRED scores. Higher is more selective (default="${quality}").
    --mapping           # 5 maps reads according to their 5' end. 3 maps reads according to their 3' end (default="${mapping}").
    --read_len_min      # Minimum length of reads included (default="${read_len_min}").
    --read_len_max      # Maximum length of reads included (default="${read_len_max}").
    --normalised        # rpm - use normalized data (reads per million raw), raw - raw counts (default="${normalised}").
    --mapped_twice      # Reads mapped once or twice (default="${mapped_twice}").
    --metagene_threshold# Raw counts around start or stop if less ignore a gene/region (default="${metagene_threshold}").
                        # If Normalized = "rpm"  MetagThreshold is readjusted by dividing normalization_factor inside program FILTER 3.
    --metagene_span     # Number of nucleotides before and after start/stop in metagene profiles (default="${metagene_span}").
    --metagene_spancorr # Number of nucleotides before and after start/stop in metagenomic profile from corrected data 
                        # (default="${metagene_span}").
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
                        # Replica is defined as order of sample names under the Group
                        # Space separated depleted samples.
    --groupB            # Space separated normal/wild type samples.

EOF
)

log "$0 $*"
# Save command line args for logging (they will be lost after utils/parse_options.sh)
run_args=$(utils/print_args.sh $0 "$@")
. utils/parse_options.sh

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

# Define a function to check if a stage is already completed
is_stage_completed() {
    local stage_number=$1
    local marker_file="stage_${stage_number}_done.marker"

    if [ -f "${marker_file}" ]; then
        log "Stage ${stage_number} has already been completed. Skipping..."
        return 0  # Indicate stage is completed
    else
        return 1  # Indicate stage is not completed
    fi
}

# Define a function to mark a stage as completed
mark_stage_completed() {
    local stage_number=$1
    local marker_file="stage_${stage_number}_done.marker"
    touch "${marker_file}"
    log "Stage ${stage_number} marked as completed."
}

if [ ${stage} -eq 1 ] && [ "${seq}" == "rna" ]; then
    stage=2 # We skip the adapter trimming stage for RNA-Seq samples
    mark_stage_completed 1
fi

if [ ${stage} -le 1 ] && [ ${stop_stage} -ge 1 ]; then
    if is_stage_completed 1; then
        echo "Remove the marker file to re-run this stage" && exit 1
    fi

    log "Step 1: Trimming adapter sequnces from the ${samples} FASTQ files in ${rawdir} using CutAdapt. 
    CutAdapt here removes the adapter sequence CTGTAGGCACCATCAAT from high-throughput sequencing reads."
    mkdir -p ${basedir}/1-Trimmed ${basedir}/1-Trimmed/Reports
    for sample in ${samples}; do
        input="${rawdir}/${sample}.fastq.gz"
        output="${basedir}/1-Trimmed/${sample}_trimmed.fastq.gz"
        cutadapt --discard-untrimmed -j ${ncpu} -a "CTGTAGGCACCATCAAT" \
            -o ${output} ${input} > ${basedir}/1-Trimmed/Reports/${sample}.txt || exit 1
    done

    # Mark stage 1 as completed
    mark_stage_completed 1
fi

if [ ${stage} -le 2 ] && [ ${stop_stage} -ge 2 ]; then
    if is_stage_completed 2; then
        echo "Remove the marker file to re-run this stage" && exit 1
    fi

    log "Step 2: Applying a quality filter on the trimmed RiboSeq or raw RNA-Seq FASTQ files. 
    Low quality reads i.e. those with a Phred score of less than 20 are discarded. Sequences outside 
    the range of ${read_len_min} and ${read_len_max} are discarded."
    ${python} local/quality_filtering.py --rawdir ${rawdir} --names "${samples}" \
        --threshold ${quality} --seq_type ${seq_type} \
        --read_len_min ${read_len_min} --read_len_max ${read_len_max}

    mark_stage_completed 2
fi

if [ ${stage} -le 3 ] && [ ${stop_stage} -ge 3 ]; then
    if is_stage_completed 3; then
        echo "Remove the marker file to re-run this stage" && exit 1
    fi

    log "Step 3: Removing non-coding mRNA sequences from the filtered FASTQ files using Bowtie2."
    mkdir -p ${basedir}/3-Subtracted ${basedir}/3-Subtracted/Reports ${basedir}/3-Subtracted/SAM
    for sample in ${samples}; do 
        ncRNA="${refdir}/Indexes/ncRNA"
        input="${basedir}/2-Filtered/${sample}_filtered.fastq"
        output="${basedir}/3-Subtracted/SAM/${sample}.sam"
        unmapped="${basedir}/3-Subtracted/${sample}.fastq"
        bowtie2 --no-unal -p ${ncpu} --un ${unmapped} -x ${ncRNA} \
            -U ${input} -S ${output} > ${basedir}/3-Subtracted/Reports/${sample}.txt || exit 1
    done
    if [ "${clean}" -eq 1 ]; then
        rm -r ${basedir}/3-Subtracted/SAM
    fi

    mark_stage_completed 3
fi

if [ ${stage} -le 4 ] && [ ${stop_stage} -ge 4 ]; then
    if is_stage_completed 4; then
        echo "Remove the marker file to re-run this stage" && exit 1
    fi

    log "Step 4: Aligning coding mRNA sequences to the S. cerevisiae reference genome using HISAT2."
    mkdir -p ${basedir}/4-Aligned ${basedir}/4-Aligned/Reports
    for sample in ${samples}; do
        genome="${refdir}/Indexes/Genome"
        input="${basedir}/3-Subtracted/${sample}.fastq"
        output="${basedir}/4-Aligned/${sample}.sam"
        hisat2 --no-unal -p ${ncpu} -k 2 --no-softclip --dta -x ${genome} \
            -U ${input} -S ${output} > ${basedir}/4-Aligned/Reports/${sample}.txt || exit 1
    done

    log "\nConverting aligned sequences from SAM to BAM with Samtools."
    for sample in ${samples}; do
        sam=${basedir}/4-Aligned/${sample}.sam
        bam=${basedir}/4-Aligned/${sample}.bam
        # convert sam to bam
        samtools view -S -b ${sam} -o ${bam}.unsorted 
        # sort converted bam
        samtools sort ${bam}.unsorted -o ${bam}
        # index sorted bam
        samtools index ${bam}
        rm ${bam}.unsorted
    done
    if [ "${clean}" -eq 1 ]; then
        rm ${basedir}/4-Aligned/*.sam
    fi

    mark_stage_completed 4
fi

if [ ${seq_type} == "rna" ]; then
    log "RNA-Seq Raw assignment and Meta-Gene tables currently not implemented. Use RiboSeq instead."
elif [ ${seq_type} == "ribo" ]; then
    if [ ${stage} -le 5 ] && [ ${stop_stage} -ge 5 ]; then
        if is_stage_completed 5; then
            echo "Remove the marker file to re-run this stage" && exit 1
        fi

        log "Step 5: Processing Raw assignments."
        ${python} local/process_raw_assignment.py --rawdir ${rawdir} --names "${samples}" \
            --mapped_twice ${mapped_twice} --mapping ${mapping} \
            --read_len_min ${read_len_min} --read_len_max ${read_len_max}

        mark_stage_completed 5
    fi

    if [ ${stage} -le 6 ] && [ ${stop_stage} -ge 6 ]; then
        if is_stage_completed 6; then
            echo "Remove the marker file to re-run this stage" && exit 1
        fi

        log "Step 6: Creating Meta-Gene Tables."
        ${python} local/create_meta_tables.py --rawdir ${rawdir} --names "${samples}" \
            --metagene_span ${metagene_span} --metagene_threshold ${metagene_threshold} \
            --mapped_twice ${mapped_twice} --mapping ${mapping} --normalised ${normalised} \
            --read_len_min ${read_len_min} --read_len_max ${read_len_max}

        mark_stage_completed 6
    fi

    if [ ${stage} -le 7 ] && [ ${stop_stage} -ge 7 ]; then
        if is_stage_completed 7; then
            echo "Remove the marker file to re-run this stage" && exit 1
        fi

        log "Step 7: Plotting normalized Meta-Gene tables."
        ${python} local/plot_meta_tables.py --rawdir ${rawdir} --names "${samples}" \
            --metagene_span ${metagene_span} --mapping ${mapping} --normalised ${normalised} \
            --read_len_min ${read_len_min} --read_len_max ${read_len_max}

        mark_stage_completed 7
    fi
fi
