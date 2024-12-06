Library amplification using SMART-9N causes repeated insertion of TSO adapter subsequences at the read ends. This tool specifically recognizes the adapter sequence from each read end and removes it, including the repeated adapter subsequences.



Compile

g++ -I ./ -std=c++17 trimTSO.cpp -lz -pthread -O3 -o trimTSO


Usage

./trimTSO -i InputF.fastq -I InputR.fastq -o output_F.fastq -O output_R.fastq -s single.fastq -a AAGCAGTGGTATCAACGCAGAGTACATGGG [-m 8 -n 1 -l 20]


gzip-hpp(https://github.com/mapbox/gzip-hpp) to (de)compress gz
