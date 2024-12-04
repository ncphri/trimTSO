Library amplification using SMART-9N causes repeated insertion of TSO adapter subsequences at the read ends. This tool specifically recognizes the adapter sequence from each read end and removes it, including the repeated adapter subsequences.

Compile

g++ -std=c++17 SMARTadapterTrimmer.cpp -lz -pthread -O3 -o SMARTadapterTrimmer

Usage

./SMARTadapterTrimmer -1 InputF.fastq -2 InputR.fastq -3 output_F.fastq -4 output_R.fastq -s single.fastq -a AAGCAGTGGTATCAACGCAGAGTACATGGG -m 8
