Library amplification using SMART-9N causes repeated insertion of TSO adapter subsequences at the read ends. This tool specifically recognizes the adapter sequence from each read end and removes it, including the repeated adapter subsequences.

This tool is also useful for adapter/primer trimming as it recognizes and trims adapter sequences from the read ends.



Compile

g++ -std=c++17 trimTSO.cpp -lz -pthread -O3 -o trimTSO


Usage

./trimTSO -i InputF.fastq(.gz) [-I InputR..fastq(.gz)] -o output_F [-O output_R -s single] -f adapters.fa [-m 8 -n 1 -l 20 -r -g]
