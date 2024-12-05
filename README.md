Library amplification using SMART-9N causes repeated insertion of TSO adapter subsequences at the read ends. This tool specifically recognizes the adapter sequence from each read end and removes it, including the repeated adapter subsequences.



Compile

g++ -I ./ -std=c++17 trimTSOgz-modified5.cpp -lz -pthread -O3 -o trimTSOgz-modified5


Usage

./trimTSO -1 InputF.fastq -2 InputR.fastq -3 output_F.fastq -4 output_R.fastq -s single.fastq -a AAGCAGTGGTATCAACGCAGAGTACATGGG -m 8


gzip-hpp(https://github.com/mapbox/gzip-hpp) to decompress gz
