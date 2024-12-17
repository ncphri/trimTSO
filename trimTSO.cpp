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
    bool recurse;  // 新しいメンバ変数

    std::mutex writeMutex;

    // バッファ用メンバ変数
    std::string out1Buffer;
    std::string out2Buffer;
    std::string singleOutBuffer;
    const size_t bufferLimit = 10 * 1024 * 1024; // 10MB

    // バッファをファイルに書き込む関数
    void flushBufferToFile(const std::string& buffer, const std::string& fileName) {
        if (!buffer.empty()) {
            std::ofstream file(fileName, std::ios_base::app);
            file << buffer;
        }
    }

    void processPair(const std::string& readname1, const std::string& read1,
                     const std::string& qual1, const std::string& readname2,
                     const std::string& read2, const std::string& qual2,
                     AdapterTrimmer& trimmer) {
        auto [trimmedSeq1, trimmedQual1] = trimmer.trimAdapters(read1, qual1);
        auto [trimmedSeq2, trimmedQual2] = trimmer.trimAdapters(read2, qual2);

        std::lock_guard<std::mutex> lock(writeMutex);

        if (trimmedSeq1.size() >= static_cast<size_t>(minReadLength) && trimmedSeq2.size() >= static_cast<size_t>(minReadLength)) {
            // 両リードがトリム可能
            out1Buffer += readname1 + "\n" + trimmedSeq1 + "\n+\n" + trimmedQual1 + "\n";
            out2Buffer += readname2 + "\n" + trimmedSeq2 + "\n+\n" + trimmedQual2 + "\n";
        } else if (trimmedSeq1.size() >= static_cast<size_t>(minReadLength)) {
            // read1のみトリム可能
            singleOutBuffer += readname1 + "\n" + trimmedSeq1 + "\n+\n" + trimmedQual1 + "\n";
        } else if (trimmedSeq2.size() >= static_cast<size_t>(minReadLength)) {
            // read2のみトリム可能
            singleOutBuffer += readname2 + "\n" + trimmedSeq2 + "\n+\n" + trimmedQual2 + "\n";
        }

        // 一定サイズを超えたらバッファをフラッシュ
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
                   bool rc = false)  // 新しいコンストラクタパラメータ 
        : inputFile1(in1), inputFile2(in2), 
          outputFile1(out1), outputFile2(out2), 
          singleOutputFile(single), 
          adapterArray(adapterSeqs), 
          matchLength(matchLen),
          minReadLength(minLen),
          maxMismatches(maxMismatchCount),
          maxMismatchCost(maxCosts),
          recurse(rc) {}  // 新しいメンバ変数の初期化

    void process() {
        AdapterTrimmer trimmer(adapterArray, matchLength, maxMismatches, maxMismatchCost, recurse);  // maxMismatchesを追加
        std::vector<std::string> readnames1, reads1, quals1;
        std::vector<std::string> readnames2, reads2, quals2;

        // Read entire contents for both input files
        readFastqFromString(inputFile1, readnames1, reads1, quals1);
        readFastqFromString(inputFile2, readnames2, reads2, quals2);

        // Ensure equal number of reads in paired files
        if (reads1.size() != reads2.size()) {
            std::cerr << "Warning: Unequal number of reads in paired input files!" << std::endl;
            return;
        }

        // Process each pair of reads
        for (size_t i = 0; i < reads1.size(); ++i) {
            processPair(readnames1[i], reads1[i], quals1[i], 
                        readnames2[i], reads2[i], quals2[i], trimmer);
        }
        // 残りのデータをファイルにフラッシュ
        flushBufferToFile(out1Buffer, outputFile1);
        flushBufferToFile(out2Buffer, outputFile2);
        flushBufferToFile(singleOutBuffer, singleOutputFile);
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

int main(int argc, char* argv[]) {
    std::string inputFile1, inputFile2, outputFile1, outputFile2, singleOutputFile, adapterFile;
    std::string decompressed_data1, decompressed_data2;
    int matchLength = 8;
    int minReadLength = 0;
    int maxMismatches = 0;
    int maxMismatchCost = 0;
    bool recurse = false;  // New flag for recursive processing
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
                          << "-i [input_forward.fastq(.gz)]\n"
                          << "-I [input_reverse.fastq(.gz)]\n"
                          << "-o [forward_output]\n"
                          << "-O [reverse_output]\n"
                          << "-s [single_output]\n"
                          << "-f [adapterfile(multiFASTA)]\n"
                          << " \n"
                          << "-------optional args-------\n"
                          << "-m min_match_length (default: 8)\n"
                          << "-l min_read_length (default: 0)\n"
                          << "-n max_mismatches_count (default: 0)\n"
                          << "-c max_mismatch_cost (This option increases computation time; default: 0)\n"
                          << "-r trim recursively (only recommended for SMART-adapter trim; default: false)\n"
                          << "-g gzip output (default: false)\n";  // Update usage
                return 1;
            default:
                std::cerr << "Unknown option: " << char(opt) << "\n";
                return 1;
        }
    }

    

    if (inputFile1.empty() || inputFile2.empty() || outputFile1.empty() || 
        outputFile2.empty() || singleOutputFile.empty() || adapterFile.empty()) {
        std::cerr << "Error: Missing required arguments." << std::endl;
        return 1;
    }
    outputFile1 = outputFile1 + ".fastq";
    outputFile2 = outputFile2 + ".fastq";
    singleOutputFile = singleOutputFile + ".fastq";

    // 確認用の出力
    std::cout << "Final Values:\n";
    std::cout << "inputFile1: " << inputFile1 << "\n";
    std::cout << "inputFile2: " << inputFile2 << "\n";
    std::cout << "outputFile1: " << outputFile1 << "\n";
    std::cout << "outputFile2: " << outputFile2 << "\n";
    std::cout << "singleOutputFile: " << singleOutputFile << "\n";
    std::cout << "adapterFile: " << adapterFile << "\n";
    std::cout << "matchLength: " << matchLength << "\n";
    std::cout << "minReadLength: " << minReadLength << "\n";
    if (maxMismatchCost > 0){
        std::cout << "maxMismatchCost: " << maxMismatchCost << "\n";
    }
    else{
        std::cout << "maxMismatches: " << maxMismatches << "\n";
    }
    std::cout << "recurse: " << (recurse ? "true" : "false") << "\n";
    std::cout << "gzipout: " << (gzipout ? "true" : "false") << "\n";


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

    const std::filesystem::path inputFiles[2] = {inputFile1, inputFile2};

for (int i = 0; i < 2; ++i) {
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
