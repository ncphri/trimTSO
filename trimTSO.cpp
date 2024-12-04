#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <getopt.h>
#include <zlib.h>

class AdapterTrimmer {
private:
    // Define TrimmedRead at the class level
    struct TrimmedRead {
        std::string sequence;
        std::string quality;
    };

    std::string adapter;
    int matchLength;

    // 逆配列の作成（相補配列）
    std::string reverseComplement(const std::string& seq) {
        std::string rc = seq;
        std::reverse(rc.begin(), rc.end());
        for (char& c : rc) {
            switch (c) {
                case 'A': c = 'T'; break;
                case 'T': c = 'A'; break;
                case 'C': c = 'G'; break;
                case 'G': c = 'C'; break;
                default: c = 'N'; break;
            }
        }
        return rc;
    }

public:
    AdapterTrimmer(const std::string& adapterSeq, int matchLen = 8)
        : adapter(adapterSeq), matchLength(matchLen) {}

    std::pair<std::string, std::string> trimAdapters(const std::string& sequence, const std::string& quality) {
        if (sequence.empty() || quality.empty()) {
            return {sequence, quality};
        }

        // 1. リードからアダプターをトリム (Forward)
        TrimmedRead forwardTrimmed = forwardTrim(sequence, quality);

        // 2. リードを相補配列に変換
        std::string rcSequence = reverseComplement(forwardTrimmed.sequence);
        std::string rcQuality = reverseString(forwardTrimmed.quality);

        // 3. 相補配列をトリム
        TrimmedRead rcTrimmed = forwardTrim(rcSequence, rcQuality);

        // 4. リードをフォワード配列に戻す（相補配列→元の配列）
        std::string finalSequence = reverseComplement(rcTrimmed.sequence);
        std::string finalQuality = reverseString(rcTrimmed.quality);

        return {finalSequence, finalQuality};
    }

private:
    // 配列の前方からのトリミング
    TrimmedRead forwardTrim(const std::string& sequence, const std::string& quality) {
        std::string trimmedSeq = sequence;
        std::string trimmedQual = quality;
        bool trimmed = true;
        while (trimmed && trimmedSeq.length() >= static_cast<size_t>(matchLength)) {
            trimmed = false;
            for (int x = matchLength; x <= static_cast<int>(adapter.length()); ++x) {
                if (trimmedSeq.length() >= static_cast<size_t>(x)) {
                    // トリミング確認（順方向アダプター）
                    std::string adapterEnd = adapter.substr(adapter.length() - x);
                    if (trimmedSeq.substr(0, x) == adapterEnd) {
                        trimmedSeq = trimmedSeq.substr(x);
                        trimmedQual = trimmedQual.substr(x);
                        trimmed = true;
                        break;
                    }
                }
            }
        }
        return {trimmedSeq, trimmedQual};
    }

    // 文字列を逆にする
    std::string reverseString(const std::string& str) {
        std::string reversed = str;
        std::reverse(reversed.begin(), reversed.end());
        return reversed;
    }
};

class FastqProcessor {
private:
    std::string inputFile1;
    std::string inputFile2;
    std::string outputFile1;
    std::string outputFile2;
    std::string singleOutputFile;
    std::string adapterSeq;
    int matchLength;
    bool gzipOutput;

    std::mutex writeMutex;

    void processPair(const std::string& readname1, const std::string& read1,
                     const std::string& qual1, const std::string& readname2,
                     const std::string& read2, const std::string& qual2,
                     AdapterTrimmer& trimmer) {
        auto [trimmedSeq1, trimmedQual1] = trimmer.trimAdapters(read1, qual1);
        auto [trimmedSeq2, trimmedQual2] = trimmer.trimAdapters(read2, qual2);

        std::lock_guard<std::mutex> lock(writeMutex);

        if (!trimmedSeq1.empty() && !trimmedSeq2.empty()) {
            // 両リードがトリム可能
            std::ofstream out1(outputFile1, std::ios_base::app);
            std::ofstream out2(outputFile2, std::ios_base::app);

            out1 << readname1 << "\n" << trimmedSeq1 << "\n+\n" << trimmedQual1 << "\n";
            out2 << readname2 << "\n" << trimmedSeq2 << "\n+\n" << trimmedQual2 << "\n";
        } else if (!trimmedSeq1.empty()) {
            // read1のみトリム可能
            std::ofstream singleOut(singleOutputFile, std::ios_base::app);
            singleOut << readname1 << "\n" << trimmedSeq1 << "\n+\n" << trimmedQual1 << "\n";
        } else if (!trimmedSeq2.empty()) {
            // read2のみトリム可能
            std::ofstream singleOut(singleOutputFile, std::ios_base::app);
            singleOut << readname2 << "\n" << trimmedSeq2 << "\n+\n" << trimmedQual2 << "\n";
        }
    }

    void readFastqFile(std::ifstream& file, std::vector<std::string>& readnames,
                       std::vector<std::string>& reads, std::vector<std::string>& quals) {
        std::string line;
        readnames.clear();
        reads.clear();
        quals.clear();

        while (readnames.size() < 10000 && std::getline(file, line)) {
            std::string sequence, plus, quality;

            // Validate @readname line
            if (line.empty() || line[0] != '@') {
                std::cerr << "Warning: Skipping malformed FASTQ entry (missing @readname)" << std::endl;
                continue;
            }
            std::string readname = line;

            if (!std::getline(file, sequence)) break;
            if (!std::getline(file, plus)) break;
            if (!std::getline(file, quality)) break;

            // Validate FASTQ format
            if (sequence.empty() || plus.empty() || quality.empty() ||
                plus[0] != '+' || sequence.length() != quality.length()) {
                std::cerr << "Warning: Skipping malformed FASTQ entry" << std::endl;
                continue;
            }

            readnames.push_back(readname);
            reads.push_back(sequence);
            quals.push_back(quality);
        }
    }

public:
    FastqProcessor(const std::string& in1, const std::string& in2,
                   const std::string& out1, const std::string& out2,
                   const std::string& single, 
                   const std::string& adapter, 
                   int matchLen = 8,
                   bool gz = false) 
        : inputFile1(in1), inputFile2(in2), 
          outputFile1(out1), outputFile2(out2), 
          singleOutputFile(single), 
          adapterSeq(adapter), 
          matchLength(matchLen),
          gzipOutput(gz) {}

    void process() {
        std::ifstream file1(inputFile1);
        std::ifstream file2(inputFile2);

        if (!file1.is_open() || !file2.is_open()) {
            std::cerr << "Error: Unable to open input files." << std::endl;
            return;
        }

        AdapterTrimmer trimmer(adapterSeq, matchLength);
        std::vector<std::string> readnames1, reads1, quals1;
        std::vector<std::string> readnames2, reads2, quals2;

        while (true) {
            reads1.clear(); reads2.clear();
            quals1.clear(); quals2.clear();
            readnames1.clear(); readnames2.clear();

            readFastqFile(file1, readnames1, reads1, quals1);
            readFastqFile(file2, readnames2, reads2, quals2);

            // Ensure equal number of reads in paired files
            if (reads1.size() != reads2.size()) {
                std::cerr << "Warning: Unequal number of reads in paired input files!" << std::endl;
                break;
            }

            if (reads1.empty() && reads2.empty()) break;

            // Process each pair of reads
            for (size_t i = 0; i < reads1.size(); ++i) {
                processPair(readnames1[i], reads1[i], quals1[i], readnames2[i], reads2[i], quals2[i], trimmer);
            }
        }
    }
};

int main(int argc, char* argv[]) {
    std::string inputFile1, inputFile2, outputFile1, outputFile2, singleOutputFile, adapterSeq;
    int matchLength = 8;
    bool gzipOutput = false;

    int opt;
    while ((opt = getopt(argc, argv, "1:2:3:4:s:a:m:g")) != -1) {
        switch (opt) {
            case '1': inputFile1 = optarg; break;
            case '2': inputFile2 = optarg; break;
            case '3': outputFile1 = optarg; break;
            case '4': outputFile2 = optarg; break;
            case 's': singleOutputFile = optarg; break;
            case 'a': adapterSeq = optarg; break;
            case 'm': matchLength = std::stoi(optarg); break;
            case 'g': gzipOutput = true; break;
            default:
                std::cerr << "Usage: " << argv[0] << " -1 input_file1 -2 input_file2 -3 output_file1 -4 output_file2 -s single_output -a adapter_seq [-m match_length]" << std::endl;
                return 1;
        }
    }

    if (inputFile1.empty() || inputFile2.empty() || outputFile1.empty() || outputFile2.empty() || singleOutputFile.empty() || adapterSeq.empty()) {
        std::cerr << "Error: Missing required arguments." << std::endl;
        return 1;
    }

    FastqProcessor processor(inputFile1, inputFile2, outputFile1, outputFile2, singleOutputFile, adapterSeq, matchLength, gzipOutput);
    processor.process();

    return 0;
}
