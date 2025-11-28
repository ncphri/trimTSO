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
#include <zlib.h>

// Structure to hold adapter information and statistics
struct AdapterInfo {
    std::string name;
    std::string sequence;
    size_t trimCount = 0;
    
    AdapterInfo(const std::string& n, const std::string& s) : name(n), sequence(s) {}
};

class AdapterTrimmer {
private:
    // Define TrimmedRead at the class level
    struct TrimmedRead {
        std::string sequence;
        std::string quality;
        int trimmedAdapterIndex = -1; // Index of adapter that was used for trimming
    };

    std::vector<AdapterInfo> adapterArray;
    int matchLength;
    int maxMismatches;
    int maxMismatchCost;
    bool recurse;
    bool partial; // New flag for partial trimming
    mutable std::mutex statsMutex; // Mutex for thread-safe statistics updates

    // 逆配列の作成（相補配列）
    std::string reverseComplement(const std::string& seq) {
        std::string rc = seq;
        std::reverse(rc.begin(), rc.end());
        for (char& c : rc) {
            c = toupper(c);
            switch (c) {
                case 'A': c = 'T'; break;
                case 'T': c = 'A'; break;
                case 'U': c = 'A'; break;
                case 'C': c = 'G'; break;
                case 'G': c = 'C'; break;
                case 'R': c = 'Y'; break; // Purine (A/G) -> Pyrimidine (T/C)
                case 'Y': c = 'R'; break; // Pyrimidine (C/T) -> Purine (G/A)
                case 'K': c = 'M'; break; // Keto (G/T) -> Amino (A/C)
                case 'M': c = 'K'; break; // Amino (A/C) -> Keto (T/G)
                case 'S': c = 'S'; break; // Strong (G/C) -> Strong (C/G)
                case 'W': c = 'W'; break; // Weak (A/T) -> Weak (T/A)
                case 'B': c = 'V'; break; // Not A (C/G/T) -> Not T (G/C/A)
                case 'D': c = 'H'; break; // Not C (A/G/T) -> Not G (T/C/A)
                case 'H': c = 'D'; break; // Not G (A/C/T) -> Not C (T/G/A)
                case 'V': c = 'B'; break; // Not T (A/C/G) -> Not A (T/G/C)
                case 'N': c = 'N'; break;
                default: c = 'N'; break;
            }
        }
        return rc;
    }

    // Check if a base matches an IUPAC code
    bool isIUPACMatch(char readBase, char adapterBase) {
        readBase = toupper(readBase);
        adapterBase = toupper(adapterBase);
        
        if (readBase == adapterBase) return true;
        if (adapterBase == 'N') return true; // Adapter N matches anything
        
        switch (adapterBase) {
            case 'R': return (readBase == 'A' || readBase == 'G');
            case 'Y': return (readBase == 'C' || readBase == 'T');
            case 'S': return (readBase == 'G' || readBase == 'C');
            case 'W': return (readBase == 'A' || readBase == 'T');
            case 'K': return (readBase == 'G' || readBase == 'T');
            case 'M': return (readBase == 'A' || readBase == 'C');
            case 'B': return (readBase != 'A'); // C, G, T
            case 'D': return (readBase != 'C'); // A, G, T
            case 'H': return (readBase != 'G'); // A, C, T
            case 'V': return (readBase != 'T'); // A, C, G
            default: return false;
        }
    }

    // 塩基間の距離を計算（他と異なる場合は1、同じ場合は0）
    int calculateBaseDistance(char base1, char base2) {
        // base1: read sequence, base2: adapter sequence (can contain IUPAC codes)
        return isIUPACMatch(base1, base2) ? 0 : 1;
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
    AdapterTrimmer(const std::vector<AdapterInfo>& adapterInfos, int matchLen = 8, int maxMismatchCount = 0, int maxCosts = 0, bool rc = false, bool p = false)
        : adapterArray(adapterInfos), matchLength(matchLen), maxMismatches(maxMismatchCount), maxMismatchCost(maxCosts), recurse(rc), partial(p) {}

    std::pair<std::string, std::string> trimAdapters(const std::string& sequence, const std::string& quality) {
        if (sequence.empty() || quality.empty()) {
            return {sequence, quality};
        }

        // 1. リードからアダプターをトリミング (Forward)
        TrimmedRead forwardTrimmed = forwardTrim(sequence, quality);

        // 2. リードを相補配列に変換
        std::string rcSequence = reverseComplement(forwardTrimmed.sequence);
        std::string rcQuality = reverseString(forwardTrimmed.quality);

        // 3. 相補配列をトリミング 
        TrimmedRead rcTrimmed = forwardTrim(rcSequence, rcQuality);

        // Update statistics based on which adapter was used
        if (forwardTrimmed.trimmedAdapterIndex != -1) {
            std::lock_guard<std::mutex> lock(statsMutex);
            adapterArray[forwardTrimmed.trimmedAdapterIndex].trimCount++;
        } else if (rcTrimmed.trimmedAdapterIndex != -1) {
            std::lock_guard<std::mutex> lock(statsMutex);
            adapterArray[rcTrimmed.trimmedAdapterIndex].trimCount++;
        }

        // 4. リードをフォワード配列に戻す（相補配列→元の配列）
        std::string finalSequence = reverseComplement(rcTrimmed.sequence);
        std::string finalQuality = reverseString(rcTrimmed.quality);

        return {finalSequence, finalQuality};
    }

    // Method to get trimming statistics
    std::vector<AdapterInfo> getTrimStats() const {
        std::lock_guard<std::mutex> lock(statsMutex);
        return adapterArray;
    }

    // Method to print trimming statistics
    void printTrimStats() const {
        std::lock_guard<std::mutex> lock(statsMutex);
        std::cout << "\n=== Adapter Trimming Statistics ===" << std::endl;
        std::cout << "Adapter Name\tSequence\tTrimmed Reads" << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        
        size_t totalTrimmed = 0;
        for (const auto& adapter : adapterArray) {
            std::cout << adapter.name << "\t" << adapter.sequence << "\t" << adapter.trimCount << std::endl;
            totalTrimmed += adapter.trimCount;
        }
        
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "Total trimmed reads: " << totalTrimmed << std::endl;
    }

    // Method to write trimming statistics to file
    void writeTrimStatsToFile(const std::string& filename) const {
        std::lock_guard<std::mutex> lock(statsMutex);
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Warning: Could not open statistics file: " << filename << std::endl;
            return;
        }
        
        outFile << "Adapter_Name\tSequence\tTrimmed_Reads\n";
        size_t totalTrimmed = 0;
        
        for (const auto& adapter : adapterArray) {
            outFile << adapter.name << "\t" << adapter.sequence << "\t" << adapter.trimCount << "\n";
            totalTrimmed += adapter.trimCount;
        }
        
        outFile << "Total\tAll_Adapters\t" << totalTrimmed << "\n";
        outFile.close();
        
        std::cout << "Trimming statistics written to: " << filename << std::endl;
    }

private:
    // 配列の前方からのトリミング（ミスマッチ許容）
    TrimmedRead forwardTrim(const std::string& sequence, const std::string& quality) {
        std::string trimmedSeq = sequence;
        std::string trimmedQual = quality;
        bool trimmed = true;
        int usedAdapterIndex = -1;
        
        for (size_t adapterIdx = 0; adapterIdx < adapterArray.size(); ++adapterIdx) {
            const auto& adapter = adapterArray[adapterIdx].sequence;
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
                        
                        bool shouldTrim = false;
                        if (maxMismatchCost > 0){
                            shouldTrim = (mismatches <= maxMismatchCost);
                        } else {
                            shouldTrim = (mismatches <= maxMismatches);
                        }
                        
                        if (shouldTrim) {
                            trimmedSeq = trimmedSeq.substr(x);
                            trimmedQual = trimmedQual.substr(x);
                            usedAdapterIndex = adapterIdx; // Record which adapter was used
                            if (recurse == true) {
                                trimmed = true;
                            }
                            break;  // 現在のadapterで一致したのでbreak
                        }
                    }
                }
            }

            // Partial adapter trimming (only if enabled)
            if (!trimmed && partial) {
                 for (int x = matchLength; x < static_cast<int>(adapter.length()); ++x) {
                    if (trimmedSeq.length() >= static_cast<size_t>(x)) {
                        // Check if the beginning of the read matches ANY part of the adapter
                        for (size_t y = 0; y <= adapter.length() - x; ++y) {
                            std::string adapterSub = adapter.substr(y, x);
                            std::string toMatch = trimmedSeq.substr(0, x);
                            
                            int mismatches = 0;
                            if (maxMismatchCost > 0){
                                mismatches = calculateEditDistance(adapterSub, toMatch);
                            }
                            else{
                                for (size_t i = 0; i < toMatch.length(); ++i) {
                                    mismatches += calculateBaseDistance(toMatch[i], adapterSub[i]);
                                }
                            }
                            
                            bool shouldTrim = false;
                            if (maxMismatchCost > 0){
                                shouldTrim = (mismatches <= maxMismatchCost);
                            } else {
                                shouldTrim = (mismatches <= maxMismatches);
                            }
                            
                            if (shouldTrim) {
                                trimmedSeq = trimmedSeq.substr(x);
                                trimmedQual = trimmedQual.substr(x);
                                usedAdapterIndex = adapterIdx;
                                if (recurse == true) {
                                    trimmed = true;
                                }
                                break; 
                            }
                        }
                        if (trimmed) break;
                    }
                 }
            }
        }
        return {trimmedSeq, trimmedQual, usedAdapterIndex};
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
    std::vector<AdapterInfo> adapterArray;
    int matchLength;
    int minReadLength;
    int maxMismatches;
    int maxMismatchCost;
    bool recurse;
    bool partial;
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

        // Helper to trim \r
        auto trim = [](std::string& s) {
            if (!s.empty() && s.back() == '\r') s.pop_back();
        };

        // Read entire file
        while (std::getline(iss, line)) {
            trim(line);
            std::string sequence, plus, quality;

            // Validate @readname line
            if (line.empty() || line[0] != '@') {
                continue;  // Skip malformed entries
            }
            std::string readname = line;

            // Check if we have enough data for a complete FASTQ entry
            if (!std::getline(iss, sequence) || 
                !std::getline(iss, plus) || 
                !std::getline(iss, quality)) {
                break;  // Stop if we can't read a complete entry
            }
            trim(sequence);
            trim(plus);
            trim(quality);

            if (plus.empty() || 
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
                   const std::vector<AdapterInfo>& adapterInfos, 
                   int matchLen = 8,
                   int minLen = 0,
                   int maxMismatchCount = 0,
                   int maxCosts = 0,
                   bool rc = false,
                   bool p = false)
        : inputFile1(in1), inputFile2(in2), 
          outputFile1(out1), outputFile2(out2), 
          singleOutputFile(single), 
          adapterArray(adapterInfos), 
          matchLength(matchLen),
          minReadLength(minLen),
          maxMismatches(maxMismatchCount),
          maxMismatchCost(maxCosts),
          recurse(rc),
          partial(p),
          isPaired(!in2.empty()) {}

    void process() {
        AdapterTrimmer trimmer(adapterArray, matchLength, maxMismatches, maxMismatchCost, recurse, partial);
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

        // Print and write statistics
        trimmer.printTrimStats();
        
        // Create statistics filename based on output file
        std::string statsFile = outputFile1.substr(0, outputFile1.find_last_of('.')) + "_trimming_stats.txt";
        trimmer.writeTrimStatsToFile(statsFile);
    }
};

// Custom wrapper functions to resolve casting issues
namespace gzip_custom {

std::string decompress(const char* data, size_t size) {
    if (!data || size == 0) return "";

    z_stream strm{};
    // 入力（全体）を最初に設定
    strm.next_in  = reinterpret_cast<Bytef*>(const_cast<char*>(data));
    strm.avail_in = static_cast<uInt>(size);

    // 16 + MAX_WBITS で gzip ヘッダ対応
    if (inflateInit2(&strm, 16 + MAX_WBITS) != Z_OK) {
        throw std::runtime_error("inflateInit2 failed");
    }

    std::string output;
    output.reserve(size * 3); // だいたいの予測（必要に応じて拡張されます）

    // 出力バッファは少し大きめ
    unsigned char outbuf[64 * 1024];

    while (true) {
        strm.next_out  = outbuf;
        strm.avail_out = static_cast<uInt>(sizeof(outbuf));

        int ret = inflate(&strm, Z_NO_FLUSH);

        // 今回の inflate 呼び出しで出た出力を追加
        size_t have = sizeof(outbuf) - strm.avail_out;
        if (have) {
            output.append(reinterpret_cast<const char*>(outbuf), have);
        }

        if (ret == Z_STREAM_END) {
            // 1 メンバーの終端に到達
            if (strm.avail_in == 0) {
                // 入力はすべて消費 -> すべてのメンバーを解凍し終えた
                break;
            }
            // まだ入力が残っている -> 次のメンバーを続けて解凍
            int r = inflateReset2(&strm, 16 + MAX_WBITS);
            if (r != Z_OK) {
                inflateEnd(&strm);
                throw std::runtime_error("inflateReset2 failed");
            }
            continue; // 次メンバーへ
        }

        if (ret == Z_OK) {
            // 継続。入力が尽きたのに Z_STREAM_END に達していない場合は壊れている可能性
            if (strm.avail_in == 0 && strm.avail_out != 0) {
                inflateEnd(&strm);
                throw std::runtime_error("Unexpected EOF in gzip stream (truncated file?)");
            }
            // まだ処理継続
            continue;
        }

        // それ以外はエラー
        inflateEnd(&strm);
        std::ostringstream oss;
        oss << "inflate failed (ret=" << ret << ")";
        throw std::runtime_error(oss.str());
    }

    inflateEnd(&strm);
    return output;
}

std::string compress(const std::string& input) {
    if (input.empty()) return "";
    z_stream strm{};
    strm.next_in  = reinterpret_cast<Bytef*>(const_cast<char*>(input.data()));
    strm.avail_in = static_cast<uInt>(input.size());
    if (deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
                     16 + MAX_WBITS, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
        throw std::runtime_error("deflateInit2 failed");
    }
    std::string output;
    unsigned char buffer[64 * 1024];
    int ret;
    do {
        strm.next_out  = buffer;
        strm.avail_out = sizeof(buffer);
        ret = deflate(&strm, strm.avail_in ? Z_NO_FLUSH : Z_FINISH);
        if (ret == Z_STREAM_ERROR) {
            deflateEnd(&strm);
            throw std::runtime_error("deflate failed");
        }
        size_t have = sizeof(buffer) - strm.avail_out;
        if (have) output.append(reinterpret_cast<const char*>(buffer), have);
    } while (ret != Z_STREAM_END);
    deflateEnd(&strm);
    return output;
}

} // namespace gzip_custom
// ---- ここまで置換 ----



// Modified function to read multi-FASTA file and extract adapter names
std::vector<AdapterInfo> readAdaptersFromFASTA(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open adapter file: " + filename);
    }

    std::vector<AdapterInfo> adapters;
    std::string line, currentSequence, currentName;
    bool inSequence = false;

    while (std::getline(file, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        if (line.empty()) continue;

        // Check if this is a header line
        if (line[0] == '>') {
            // If we were previously building a sequence, add it to adapters
            if (!currentSequence.empty() && !currentName.empty()) {
                // Remove any whitespace from the sequence
                currentSequence.erase(
                    std::remove_if(currentSequence.begin(), currentSequence.end(), ::isspace), 
                    currentSequence.end()
                );
                adapters.emplace_back(currentName, currentSequence);
                currentSequence.clear();
            }
            
            // Extract adapter name (remove '>' and take everything up to first space)
            currentName = line.substr(1); // Remove '>'
            size_t spacePos = currentName.find_first_of(" \t");
            if (spacePos != std::string::npos) {
                currentName = currentName.substr(0, spacePos);
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
    if (!currentSequence.empty() && !currentName.empty()) {
        // Remove any whitespace from the sequence
        currentSequence.erase(
            std::remove_if(currentSequence.begin(), currentSequence.end(), ::isspace), 
            currentSequence.end()
        );
        adapters.emplace_back(currentName, currentSequence);
    }

    return adapters;
}

bool is_gzip_compressed(const char* data, size_t size) {
    // gzipのマジックナンバーは 0x1f 0x8b
    return size >= 2 && static_cast<unsigned char>(data[0]) == 0x1f && static_cast<unsigned char>(data[1]) == 0x8b;
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
    bool partial = false;
    bool gzipout = false;
    std::vector<AdapterInfo> adapterArray;

    int opt;
    while ((opt = getopt(argc, argv, "i:I:o:O:s:f:m:l:n:c:rpgh")) != -1) {
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
            case 'p':
                partial = true;
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
                          << "-p trim partial adapters (default: false)\n"
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
              << "Partial: " << (partial ? "true" : "false") << "\n"
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
        std::cout << adapter.name << ": " << adapter.sequence << std::endl;
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
    if (is_gzip_compressed(buffer.data(), buffer.size())) {
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
                             maxMismatches, maxMismatchCost, recurse, partial);  // Add maxMismatches to constructor
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
        if (!outputFile1.empty() && std::filesystem::exists(outputFile1)) {
            compressAndWrite(outputFile1);
        }

        // outputFile2の圧縮
        if (!outputFile2.empty() && std::filesystem::exists(outputFile2)) {
            compressAndWrite(outputFile2);
        }

        // singleOutputFileの圧縮
        if (!singleOutputFile.empty() && std::filesystem::exists(singleOutputFile)) {
            compressAndWrite(singleOutputFile);
        }

    }

    return 0;
}
