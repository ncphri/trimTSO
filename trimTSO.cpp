#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <algorithm>
#include <mutex>
#include <getopt.h>
#include <zlib.h>
#include <unistd.h>
#include <gzip/compress.hpp>
#include <gzip/config.hpp>
#include <gzip/decompress.hpp>
#include <gzip/utils.hpp>
#include <gzip/version.hpp>

class AdapterTrimmer {
private:
    // Define TrimmedRead at the class level
    struct TrimmedRead {
        std::string sequence;
        std::string quality;
    };

    std::vector<std::string> adapterArray;
    int matchLength;
    int maxMismatches;
    int maxMismatchCost;
    bool recurse;

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

    // 塩基間の距離を計算（他と異なる場合は1、同じ場合は0）
    int calculateBaseDistance(char base1, char base2) {
        // 塩基が異なる場合は1、同じ場合は0を返す
        return (base1 != base2) ? 1 : 0;
    }

    int calculateEditDistance(const std::string& seq1, const std::string& seq2) {
        size_t len1 = seq1.length();
        size_t len2 = seq2.length();
    
        // DPテーブルを初期化
        std::vector<std::vector<int>> dp(len1 + 1, std::vector<int>(len2 + 1, 0));
    
        // 初期条件：空文字列への変換コスト
        for (size_t i = 0; i <= len1; ++i) dp[i][0] = i; // デリーション
        for (size_t j = 0; j <= len2; ++j) dp[0][j] = j; // インサーション
    
        // DPテーブルの更新
        for (size_t i = 1; i <= len1; ++i) {
            for (size_t j = 1; j <= len2; ++j) {
                int cost = calculateBaseDistance(seq1[i - 1], seq2[j - 1]); // ミスマッチコスト
                dp[i][j] = std::min({ 
                    dp[i - 1][j] + 1,     // デリーション（seq1の1文字を削除）
                    dp[i][j - 1] + 1,      // インサーション（seq2に1文字挿入）
                    dp[i - 1][j - 1] + cost // 置換（ミスマッチなら1, 一致なら0）
                });
            }
        }
    
        return dp[len1][len2]; // 編集距離を返す
    }

public:
    AdapterTrimmer(const std::vector<std::string>& adapterSeqs, int matchLen = 8, int maxMismatchCount = 0, int maxCosts = 0, bool rc = false)
        : adapterArray(adapterSeqs), matchLength(matchLen), maxMismatches(maxMismatchCount), maxMismatchCost(maxCosts), recurse(rc) {}

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
    // 配列の前方からのトリミング（ミスマッチ許容）
    TrimmedRead forwardTrim(const std::string& sequence, const std::string& quality) {
        std::string trimmedSeq = sequence;
        std::string trimmedQual = quality;
        bool trimmed = true;
        for (const auto& adapter : adapterArray) {
    trimmed = true;  // adapterごとにトリミングを再度有効にする
    while (trimmed && trimmedSeq.length() >= static_cast<size_t>(matchLength)) {
        trimmed = false;  // ここで false にリセットする
        for (int x = static_cast<int>(adapter.length()); x >= matchLength; --x) {
            if (trimmedSeq.length() >= static_cast<size_t>(x)) {
                // トリミング確認（順方向アダプター）
                std::string adapterEnd = adapter.substr(adapter.length() - x);
                std::string toMatch = trimmedSeq.substr(0, x);

                // ミスマッチチェック
                int mismatches = 0;
                if (maxMismatchCost > 0){
                    mismatches = calculateEditDistance(adapterEnd, toMatch);
                }
                else{
                    for (size_t i = 0; i < toMatch.length(); ++i) {
                    mismatches += calculateBaseDistance(toMatch[i], adapterEnd[i]);
                    }
                }
                
                if (maxMismatchCost > 0){
                    if (mismatches <= maxMismatchCost) {
                        trimmedSeq = trimmedSeq.substr(x);
                        trimmedQual = trimmedQual.substr(x);
                        if (recurse == true) {
                            trimmed = true;
                        }
                        break;  // 現在のadapterで一致したのでbreak
                    }
                }
                else{
                    if (mismatches <= maxMismatches) {
                        trimmedSeq = trimmedSeq.substr(x);
                        trimmedQual = trimmedQual.substr(x);
                        if (recurse == true) {
                            trimmed = true;
                        }
                        break;  // 現在のadapterで一致したのでbreak
                    }
                }
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
    std::vector<std::string> adapterArray;
    int matchLength;
    int minReadLength;
    int maxMismatches;
    int maxMismatchCost;
    bool recurse;
    bool isPaired;  // New flag to indicate if processing paired-end data

    std::mutex writeMutex;

    // Buffer variables
    std::string out1Buffer;
    std::string out2Buffer;
    std::string singleOutBuffer;
    const size_t bufferLimit = 10 * 1024 * 1024; // 10MB

    void flushBufferToFile(const std::string& buffer, const std::string& fileName) {
        if (!buffer.empty()) {
            std::ofstream file(fileName, std::ios_base::app);
            file << buffer;
        }
    }

    // New method for processing single reads
    void processSingleRead(const std::string& readname, const std::string& read,
                          const std::string& qual, AdapterTrimmer& trimmer) {
        auto [trimmedSeq, trimmedQual] = trimmer.trimAdapters(read, qual);

        std::lock_guard<std::mutex> lock(writeMutex);

        if (trimmedSeq.size() >= static_cast<size_t>(minReadLength)) {
            out1Buffer += readname + "\n" + trimmedSeq + "\n+\n" + trimmedQual + "\n";

            if (out1Buffer.size() >= bufferLimit) {
                flushBufferToFile(out1Buffer, outputFile1);
                out1Buffer.clear();
            }
        }
    }

    void processPair(const std::string& readname1, const std::string& read1,
                     const std::string& qual1, const std::string& readname2,
                     const std::string& read2, const std::string& qual2,
                     AdapterTrimmer& trimmer) {
        auto [trimmedSeq1, trimmedQual1] = trimmer.trimAdapters(read1, qual1);
        auto [trimmedSeq2, trimmedQual2] = trimmer.trimAdapters(read2, qual2);

        std::lock_guard<std::mutex> lock(writeMutex);

        if (trimmedSeq1.size() >= static_cast<size_t>(minReadLength) && 
            trimmedSeq2.size() >= static_cast<size_t>(minReadLength)) {
            out1Buffer += readname1 + "\n" + trimmedSeq1 + "\n+\n" + trimmedQual1 + "\n";
            out2Buffer += readname2 + "\n" + trimmedSeq2 + "\n+\n" + trimmedQual2 + "\n";
        } else if (trimmedSeq1.size() >= static_cast<size_t>(minReadLength)) {
            singleOutBuffer += readname1 + "\n" + trimmedSeq1 + "\n+\n" + trimmedQual1 + "\n";
        } else if (trimmedSeq2.size() >= static_cast<size_t>(minReadLength)) {
            singleOutBuffer += readname2 + "\n" + trimmedSeq2 + "\n+\n" + trimmedQual2 + "\n";
        }

        // Buffer management
        if (out1Buffer.size() >= bufferLimit) {
            flushBufferToFile(out1Buffer, outputFile1);
            out1Buffer.clear();
        }
        if (out2Buffer.size() >= bufferLimit) {
            flushBufferToFile(out2Buffer, outputFile2);
            out2Buffer.clear();
        }
        if (singleOutBuffer.size() >= bufferLimit) {
            flushBufferToFile(singleOutBuffer, singleOutputFile);
            singleOutBuffer.clear();
        }
    }

    void readFastqFromString(const std::string& fastqContent, 
                       std::vector<std::string>& readnames,
                       std::vector<std::string>& reads, 
                       std::vector<std::string>& quals) {
        std::istringstream iss(fastqContent);
        std::string line;
        readnames.clear();
        reads.clear();
        quals.clear();

        // Read entire file
        while (std::getline(iss, line)) {
            std::string sequence, plus, quality;

            // Validate @readname line
            if (line.empty() || line[0] != '@') {
                continue;  // Skip malformed entries
            }
            std::string readname = line;

            // Check if we have enough data for a complete FASTQ entry
            if (!std::getline(iss, sequence) || 
                !std::getline(iss, plus) || 
                !std::getline(iss, quality) ||
                plus.empty() || 
                plus[0] != '+' || 
                sequence.length() != quality.length()) {
                break;  // Stop if we can't read a complete entry
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
                   const std::vector<std::string>& adapterSeqs, 
                   int matchLen = 8,
                   int minLen = 0,
                   int maxMismatchCount = 0,
                   int maxCosts = 0,
                   bool rc = false)
        : inputFile1(in1), inputFile2(in2), 
          outputFile1(out1), outputFile2(out2), 
          singleOutputFile(single), 
          adapterArray(adapterSeqs), 
          matchLength(matchLen),
          minReadLength(minLen),
          maxMismatches(maxMismatchCount),
          maxMismatchCost(maxCosts),
          recurse(rc),
          isPaired(!in2.empty()) {}

    void process() {
        AdapterTrimmer trimmer(adapterArray, matchLength, maxMismatches, maxMismatchCost, recurse);
        std::vector<std::string> readnames1, reads1, quals1;
        std::vector<std::string> readnames2, reads2, quals2;

        // Read first input file
        readFastqFromString(inputFile1, readnames1, reads1, quals1);

        if (isPaired) {
            // Process paired-end data
            readFastqFromString(inputFile2, readnames2, reads2, quals2);

            if (reads1.size() != reads2.size()) {
                std::cerr << "Warning: Unequal number of reads in paired input files!" << std::endl;
                return;
            }

            for (size_t i = 0; i < reads1.size(); ++i) {
                processPair(readnames1[i], reads1[i], quals1[i], 
                           readnames2[i], reads2[i], quals2[i], trimmer);
            }

            // Flush remaining paired-end data
            flushBufferToFile(out1Buffer, outputFile1);
            flushBufferToFile(out2Buffer, outputFile2);
            flushBufferToFile(singleOutBuffer, singleOutputFile);
        } else {
            // Process single-end data
            for (size_t i = 0; i < reads1.size(); ++i) {
                processSingleRead(readnames1[i], reads1[i], quals1[i], trimmer);
            }

            // Flush remaining single-end data
            flushBufferToFile(out1Buffer, outputFile1);
        }
    }
};

// Custom wrapper functions to resolve casting issues
namespace gzip_custom {
    std::string decompress(const char* data, size_t size) {
        std::string output;
        gzip::Decompressor decomp;
        
        // Create a non-const copy to use with reinterpret_cast
        std::vector<char> data_copy(data, data + size);
        
        decomp.decompress(output, data_copy.data(), data_copy.size());
        return output;
    }

    std::string compress(const std::string& input) {
        std::string output;
        gzip::Compressor comp;
        
        // Create a non-const copy to use with reinterpret_cast
        std::string input_copy = input;
        
        comp.compress(output, input_copy.data(), input_copy.size());
        return output;
    }
}

// Custom function to read multi-FASTA file
std::vector<std::string> readAdaptersFromFASTA(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open adapter file: " + filename);
    }

    std::vector<std::string> adapters;
    std::string line, currentSequence;
    bool inSequence = false;

    while (std::getline(file, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        if (line.empty()) continue;

        // Check if this is a header line
        if (line[0] == '>') {
            // If we were previously building a sequence, add it to adapters
            if (!currentSequence.empty()) {
                // Remove any whitespace from the sequence
                currentSequence.erase(
                    std::remove_if(currentSequence.begin(), currentSequence.end(), ::isspace), 
                    currentSequence.end()
                );
                adapters.push_back(currentSequence);
                currentSequence.clear();
            }
            inSequence = true;
            continue;
        }

        // If we're in a sequence section, accumulate the sequence
        if (inSequence) {
            currentSequence += line;
        }
    }

    // Add the last sequence if not empty
    if (!currentSequence.empty()) {
        // Remove any whitespace from the sequence
        currentSequence.erase(
            std::remove_if(currentSequence.begin(), currentSequence.end(), ::isspace), 
            currentSequence.end()
        );
        adapters.push_back(currentSequence);
    }

    return adapters;
}

// Main function modifications
int main(int argc, char* argv[]) {
    std::string inputFile1, inputFile2, outputFile1, outputFile2, singleOutputFile, adapterFile;
    std::string decompressed_data1, decompressed_data2;
    int matchLength = 8;
    int minReadLength = 0;
    int maxMismatches = 0;
    int maxMismatchCost = 0;
    bool recurse = false;
    bool gzipout = false;
    std::vector<std::string> adapterArray;

    int opt;
    while ((opt = getopt(argc, argv, "i:I:o:O:s:f:m:l:n:c:rgh")) != -1) {
        switch (opt) {
            case 'i':
                inputFile1 = optarg;
                break;
            case 'I':
                inputFile2 = optarg;
                break;
            case 'o':
                outputFile1 = optarg;
                break;
            case 'O':
                outputFile2 = optarg;
                break;
            case 's':
                singleOutputFile = optarg;
                break;
            case 'f':
                adapterFile = optarg;
                break;
            case 'm':
                matchLength = std::stoi(optarg);
                break;
            case 'l':
                minReadLength = std::stoi(optarg);
                break;
            case 'n':
                maxMismatches = std::stoi(optarg);
                break;
            case 'c':
                maxMismatchCost = std::stoi(optarg);
                break;
            case 'r':
                recurse = true;
                break;
            case 'g':
                gzipout = true;
                break;
            case 'h':
                std::cout << "Usage: ./trimTSO [args]\n"
                          << "-------required args-------\n"
                          << "-i [input.fastq(.gz)]\n"
                          << "-o [output]\n"
                          << "-f [adapterfile(multiFASTA)]\n"
                          << "\n"
                          << "-------optional args-------\n"
                          << "-I [input_reverse.fastq(.gz)] (for paired-end data)\n"
                          << "-O [reverse_output] (required if -I is specified)\n"
                          << "-s [single_output] (required for paired-end data)\n"
                          << "-m min_match_length (default: 8)\n"
                          << "-l min_read_length (default: 0)\n"
                          << "-n max_mismatches_count (default: 0)\n"
                          << "-c max_mismatch_cost (increases computation time; default: 0)\n"
                          << "-r trim recursively (only for SMART-adapter trim; default: false)\n"
                          << "-g gzip output (default: false)\n";
                return 1;
            default:
                std::cerr << "Unknown option: " << char(opt) << "\n";
                return 1;
        }
    }

    // Validate required arguments
    if (inputFile1.empty() || outputFile1.empty() || adapterFile.empty()) {
        std::cerr << "Error: Missing required arguments (-i, -o, -f)" << std::endl;
        return 1;
    }

    // Validate paired-end arguments
    if (!inputFile2.empty() && outputFile2.empty()) {
        std::cerr << "Error: -O (reverse output) is required when -I is specified" << std::endl;
        return 1;
    }

    if (!inputFile2.empty() && singleOutputFile.empty()) {
        std::cerr << "Error: -s (single output) is required for paired-end data" << std::endl;
        return 1;
    }

    // Add .fastq extension to output files
    outputFile1 = outputFile1 + ".fastq";
    if (!outputFile2.empty()) {
        outputFile2 = outputFile2 + ".fastq";
        singleOutputFile = singleOutputFile + ".fastq";
    }

    // Print configuration
    std::cout << "Configuration:\n"
              << "Input file 1: " << inputFile1 << "\n";
    if (!inputFile2.empty()) {
        std::cout << "Input file 2: " << inputFile2 << "\n"
                  << "Output file 2: " << outputFile2 << "\n"
                  << "Single output: " << singleOutputFile << "\n";
    }
    std::cout << "Output file 1: " << outputFile1 << "\n"
              << "Adapter file: " << adapterFile << "\n"
              << "Match length: " << matchLength << "\n"
              << "Min read length: " << minReadLength << "\n";
    if (maxMismatchCost > 0) {
        std::cout << "Max mismatch cost: " << maxMismatchCost << "\n";
    } else {
        std::cout << "Max mismatches: " << maxMismatches << "\n";
    }
    std::cout << "Recursive: " << (recurse ? "true" : "false") << "\n"
              << "Gzip output: " << (gzipout ? "true" : "false") << "\n";

    // Read adapter sequences from multi-FASTA file
    try {
        adapterArray = readAdaptersFromFASTA(adapterFile);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // Verify at least one adapter was read
    if (adapterArray.empty()) {
        std::cerr << "Error: No adapter sequences found in the input file." << std::endl;
        return 1;
    }

    // Print out the adapter sequences (for debugging/verification)
    std::cout << "Loaded " << adapterArray.size() << " adapter sequences:" << std::endl;
    for (const auto& adapter : adapterArray) {
        std::cout << adapter << std::endl;
    }

// 入力ファイルを処理する部分を修正
int numFiles = inputFile2.empty() ? 1 : 2;  // 処理するファイル数を決定
const std::filesystem::path inputFiles[2] = {inputFile1, inputFile2};

for (int i = 0; i < numFiles; ++i) {
    size_t size = std::filesystem::file_size(inputFiles[i]);
    std::vector<char> buffer(size);

    // ファイルをバイナリモードで読み込む
    std::ifstream ifs(inputFiles[i], std::ios_base::binary);
    if (!ifs) {
        std::cerr << "Failed to open file: " << inputFiles[i] << std::endl;
        return 1;
    }
    ifs.read(buffer.data(), size);
    ifs.close();

    // 解凍結果を格納する変数
    std::string& decompressed_data = (i == 0) ? decompressed_data1 : decompressed_data2;

    // gzip圧縮されているか判定
    if (gzip::is_compressed(buffer.data(), buffer.size())) {
        std::cout << "File " << inputFiles[i] << " is gzip-compressed. Decompressing..." << std::endl;
        decompressed_data = gzip_custom::decompress(buffer.data(), buffer.size());
        std::cout << "decompressed" << std::endl;
    } else {
        std::cout << "File " << inputFiles[i] << " is not gzip-compressed. Reading as is..." << std::endl;
        decompressed_data.assign(buffer.begin(), buffer.end());  // そのまま読み込み
        std::cout << "done" << std::endl;
    }
}

    // Modify the processor creation to pass maxMismatches
    std::cout << "Processing..." << std::endl;
    FastqProcessor processor(decompressed_data1, decompressed_data2, 
                             outputFile1, outputFile2, singleOutputFile, 
                             adapterArray, matchLength, minReadLength, 
                             maxMismatches, maxMismatchCost, recurse);  // Add maxMismatches to constructor
    processor.process();
    std::cout << "Completed" << std::endl;

    // 圧縮処理を行う関数
    auto compressAndWrite = [](const std::string& outputFile) {
        size_t size = std::filesystem::file_size(outputFile);
        std::vector<char> buffer(size);

        std::ifstream ifs(outputFile, std::ios_base::binary);
        ifs.read(buffer.data(), size);
        ifs.close();

        std::string compressed_data = gzip_custom::compress(std::string(buffer.data(), size));

        std::ofstream ofs(outputFile + ".gz", std::ios_base::binary);
        ofs << compressed_data;
        ofs.close();

        std::filesystem::remove(outputFile);
    };
    
    if (gzipout == true) {
        std::cout << "Compressing output files" << std::endl;
        // outputFile1の圧縮
        compressAndWrite(outputFile1);

        // outputFile2の圧縮
        compressAndWrite(outputFile2);

        // singleOutputFileの圧縮
        compressAndWrite(singleOutputFile);
    }

    return 0;
}
