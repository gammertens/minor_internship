#!/bin/bash

process_bam() {
    file="$1"
    samtools view "$file" -h | gawk '/^@/ || (/UB:/ && /CB:/)' | samtools view -h -b > "${file%.*}_processed.bam"
}

export -f process_bam

find . -name '*.bam' -print0 | parallel -0 -j6 process_bam {}