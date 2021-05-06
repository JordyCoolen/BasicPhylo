#!/bin/bash

# USAGE
# bash run_batch_docker.sh <inputpath> <file_extension> <threads> <image>
# Example:
# bash run_batch_docker.sh /workflow/input/ _001.fastq.gz 8 jonovox/easyseq_covid19:latest
# <file_extension> most common _001.fastq.gz

for fname in ${1}/*_R1${2}
do
    base=${fname##*/}
    base=${base%_R1*}
    echo "${base}_R1${2}"
    echo "${base}_R2${2}"

    docker run -it --rm --mount type=bind,source=${PWD},target=/workflow \
    --mount type=bind,source=${1},target=/workflow/input \
    ${4} nextflow run BasicPhylo.nf \
    --reads "/workflow/input/${base}_R{1,2}${2}" --outDir /workflow/output/ \
    --threads ${3} --bootstrap 1000 --proportion 0.3 -resume
done