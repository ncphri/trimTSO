#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <mutex>
#include <getopt.h>
#include <zlib.h>
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

    std::string adapter;
    int matchLength;
    int maxMismatches;

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
    int calculateDistance(char a, char b) {
        return (a != b) ? 1 : 0;
    }

public:
    AdapterTrimmer(const std::string& adapterSeq, int matchLen = 8, int maxMismatchCount = 0)
        : adapter(adapterSeq), matchLength(matchLen), maxMismatches(maxMismatchCount) {}

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
        while (trimmed && trimmedSeq.length() >= static_cast<size_t>(matchLength)) {
            trimmed = false;
            for (int x = matchLength; x <= static_cast<int>(adapter.length()); ++x) {
                if (trimmedSeq.length() >= static_cast<size_t>(x)) {
                    // トリミング確認（順方向アダプター）
                    std::string adapterEnd = adapter.substr(adapter.length() - x);
                    std::string toMatch = trimmedSeq.substr(0, x);

                    // ミスマッチチェック
                    int mismatches = 0;
                    for (size_t i = 0; i < toMatch.length(); ++i) {
                        mismatches += calculateDistance(toMatch[i], adapterEnd[i]);
                    }

                    if (mismatches <= maxMismatches) {
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
    int minReadLength;
    int maxMismatches;  // 新しいメンバ変数

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
                   const std::string& adapter, 
                   int matchLen = 8,
                   int minLen = 0,
                   int maxMismatchCount = 0)  // 新しいコンストラクタパラメータ 
        : inputFile1(in1), inputFile2(in2), 
          outputFile1(out1), outputFile2(out2), 
          singleOutputFile(single), 
          adapterSeq(adapter), 
          matchLength(matchLen),
          minReadLength(minLen),
          maxMismatches(maxMismatchCount) {}  // 新しいメンバ変数の初期化

    void process() {
        AdapterTrimmer trimmer(adapterSeq, matchLength, maxMismatches);  // maxMismatchesを追加
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

int main(int argc, char* argv[]) {
    std::string inputFile1, inputFile2, outputFile1, outputFile2, singleOutputFile, adapterSeq;
    std::string decompressed_data1, decompressed_data2;
    int matchLength = 8;
    int minReadLength = 0;
    int maxMismatches = 0;  // New parameter for mismatch tolerance

    int opt;
    while ((opt = getopt(argc, argv, "i:I:o:O:s:a:m:l:n:h")) != -1) {
        switch (opt) {
            case 'i': inputFile1 = optarg; break;
            case 'I': inputFile2 = optarg; break;
            case 'o': outputFile1 = optarg; break;
            case 'O': outputFile2 = optarg; break;
            case 's': singleOutputFile = optarg; break;
            case 'a': adapterSeq = optarg; break;
            case 'm': matchLength = std::stoi(optarg); break;
            case 'l': minReadLength = std::stoi(optarg); break;
            case 'n': maxMismatches = std::stoi(optarg); break;  // New option for mismatch tolerance
            case 'h':
                std::cout << "Usage: ./trimTSO [args]\n"
                          << "-------required args-------\n"
                          << "-i [input_forward.fastq.gz]\n"
                          << "-I [input_reverse.fastq.gz]\n"
                          << "-o [forward_output]\n"
                          << "-O [reverse_output]\n"
                          << "-s [single_output]\n"
                          << "-a [adapter_seq]\n"
                          << " \n"
                          << "-------optional args-------\n"
                          << "-m min_match_length (default: 8)\n"
                          << "-l min_read_length (default: 0)\n"
                          << "-n max_mismatches (default: 0)\n";  // Update usage
                return 1;
            default:std::cout << "type -h for usage\n";
                return 1;
                
        }
    }

    if (inputFile1.empty() || inputFile2.empty() || outputFile1.empty() || outputFile2.empty() || singleOutputFile.empty() || adapterSeq.empty()) {
        std::cerr << "Error: Missing required arguments." << std::endl;
        return 1;
    }

    const std::filesystem::path inputFiles[2] = {inputFile1, inputFile2};

    for (int i = 0; i < 2; ++i) {
        size_t size = std::filesystem::file_size(inputFiles[i]);
        std::vector<char> buffer(size);

        std::ifstream ifs(inputFiles[i], std::ios_base::binary);
        ifs.read(buffer.data(), size);
        ifs.close();

        if (!gzip::is_compressed(buffer.data(), buffer.size())) {
            std::cerr << "File " << inputFiles[i] << " is not in gzip format" << std::endl;
            return 1;
        }

        std::string& decompressed_data = (i == 0) ? decompressed_data1 : decompressed_data2;
        decompressed_data = gzip_custom::decompress(buffer.data(), buffer.size());
    }

    // Modify the processor creation to pass maxMismatches
    FastqProcessor processor(decompressed_data1, decompressed_data2, 
                             outputFile1, outputFile2, singleOutputFile, 
                             adapterSeq, matchLength, minReadLength, 
                             maxMismatches);  // Add maxMismatches to constructor
    processor.process();

    // 圧縮処理を行う関数
    auto compressAndWrite = [](const std::string& outputFile) {
        size_t size = std::filesystem::file_size(outputFile);
        std::vector<char> buffer(size);

        std::ifstream ifs(outputFile, std::ios_base::binary);
        ifs.read(buffer.data(), size);
        ifs.close();

        std::string compressed_data = gzip_custom::compress(std::string(buffer.data(), size));

        std::ofstream ofs(outputFile + ".fastq.gz", std::ios_base::binary);
        ofs << compressed_data;
        ofs.close();

        std::filesystem::remove(outputFile);
    };

    // outputFile1の圧縮
    compressAndWrite(outputFile1);

    // outputFile2の圧縮
    compressAndWrite(outputFile2);

    // singleOutputFileの圧縮
    compressAndWrite(singleOutputFile);

    return 0;
}
