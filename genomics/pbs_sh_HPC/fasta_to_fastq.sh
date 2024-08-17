#!/bin/bash

# Loop
# https://nf-co.re/modules/art_illumina
# https://manpages.debian.org/testing/art-nextgen-simulation-tools/art_illumina.1.en.html

for fasta_file in *.fasta; do
    base_name=$(basename "$fasta_file" .fasta)

    forward_fastq="${base_name}_1.fq"
    reverse_fastq="${base_name}_2.fq"

    art_illumina -i "$fasta_file" -p -l 500 -f 20 -m 500 -s 10 -o "inputs_fastq_ariba/${base_name}_"
    
    mv "${base_name}_1.fq" "$forward_fastq"
    mv "$${base_name}_2.fq" "$reverse_fastq"
    
    echo "Processed $fasta_file -> $forward_fastq, $reverse_fastq"
done


