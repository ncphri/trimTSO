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
#include <thread>
#include <getopt.h>
#include <zlib.h>
#include <unistd.h>
#include <zlib.h>
#include <cstring>

// Structure to hold adapter information and statistics
struct AdapterInfo {
    std::string name;
    std::string sequence;
    size_t trimCount = 0;
    
    AdapterInfo(const std::string& n, const std::string& s) : name(n), sequence(s) {}
};

class AdapterTrimmer {
private:
    struct TrimmedRead {
        std::string sequence;
        std::string quality;
        int trimmedAdapterIndex = -1;
    };

    std::vector<AdapterInfo> adapterArray;
    int matchLength;
    int maxMismatches;
    int maxMismatchCost;
    bool recurse;
    bool partial;
    mutable std::mutex statsMutex;

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
                case 'R': c = 'Y'; break;
                case 'Y': c = 'R'; break;
                case 'K': c = 'M'; break;
                case 'M': c = 'K'; break;
                case 'S': c = 'S'; break;
                case 'W': c = 'W'; break;
                case 'B': c = 'V'; break;
                case 'D': c = 'H'; break;
                case 'H': c = 'D'; break;
                case 'V': c = 'B'; break;
                case 'N': c = 'N'; break;
                default: c = 'N'; break;
            }
        }
        return rc;
    }

    bool isIUPACMatch(char readBase, char adapterBase) {
        readBase = toupper(readBase);
        adapterBase = toupper(adapterBase);
        if (readBase == adapterBase) return true;
        if (adapterBase == 'N') return true;
        switch (adapterBase) {
            case 'R': return (readBase == 'A' || readBase == 'G');
            case 'Y': return (readBase == 'C' || readBase == 'T');
            case 'S': return (readBase == 'G' || readBase == 'C');
            case 'W': return (readBase == 'A' || readBase == 'T');
            case 'K': return (readBase == 'G' || readBase == 'T');
            case 'M': return (readBase == 'A' || readBase == 'C');
            case 'B': return (readBase != 'A');
            case 'D': return (readBase != 'C');
            case 'H': return (readBase != 'G');
            case 'V': return (readBase != 'T');
            default: return false;
        }
    }

    int calculateBaseDistance(char base1, char base2) {
        return isIUPACMatch(base1, base2) ? 0 : 1;
    }

    int calculateEditDistance(const std::string& seq1, const std::string& seq2) {
        size_t len1 = seq1.length();
        size_t len2 = seq2.length();
        std::vector<std::vector<int>> dp(len1 + 1, std::vector<int>(len2 + 1, 0));
        for (size_t i = 0; i <= len1; ++i) dp[i][0] = i;
        for (size_t j = 0; j <= len2; ++j) dp[0][j] = j;
        for (size_t i = 1; i <= len1; ++i) {
            for (size_t j = 1; j <= len2; ++j) {
                int cost = calculateBaseDistance(seq2[j - 1], seq1[i - 1]);
                dp[i][j] = std::min({
                    dp[i - 1][j] + 1,
                    dp[i][j - 1] + 1,
                    dp[i - 1][j - 1] + cost
                });
            }
        }
        return dp[len1][len2];
    }

public:
    AdapterTrimmer(const std::vector<AdapterInfo>& adapterInfos, int matchLen = 8, int maxMismatchCount = 0, int maxCosts = 0, bool rc = false, bool p = false)
        : adapterArray(adapterInfos), matchLength(matchLen), maxMismatches(maxMismatchCount), maxMismatchCost(maxCosts), recurse(rc), partial(p) {}

    std::pair<std::string, std::string> trimAdapters(const std::string& sequence, const std::string& quality) {
        if (sequence.empty() || quality.empty()) return {sequence, quality};
        TrimmedRead forwardTrimmed = forwardTrim(sequence, quality);
        std::string rcSequence = reverseComplement(forwardTrimmed.sequence);
        std::string rcQuality = reverseString(forwardTrimmed.quality);
        TrimmedRead rcTrimmed = forwardTrim(rcSequence, rcQuality);
        if (forwardTrimmed.trimmedAdapterIndex != -1) {
            std::lock_guard<std::mutex> lock(statsMutex);
            adapterArray[forwardTrimmed.trimmedAdapterIndex].trimCount++;
        } else if (rcTrimmed.trimmedAdapterIndex != -1) {
            std::lock_guard<std::mutex> lock(statsMutex);
            adapterArray[rcTrimmed.trimmedAdapterIndex].trimCount++;
        }
        std::string finalSequence = reverseComplement(rcTrimmed.sequence);
        std::string finalQuality = reverseString(rcTrimmed.quality);
        return {finalSequence, finalQuality};
    }

    std::vector<AdapterInfo> getTrimStats() const {
        std::lock_guard<std::mutex> lock(statsMutex);
        return adapterArray;
    }

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
    TrimmedRead forwardTrim(const std::string& sequence, const std::string& quality) {
        std::string trimmedSeq = sequence;
        std::string trimmedQual = quality;
        bool trimmed = true;
        int usedAdapterIndex = -1;
        for (size_t adapterIdx = 0; adapterIdx < adapterArray.size(); ++adapterIdx) {
            const auto& adapter = adapterArray[adapterIdx].sequence;
            trimmed = true;
            while (trimmed && trimmedSeq.length() >= static_cast<size_t>(matchLength)) {
                trimmed = false;
                for (int x = static_cast<int>(adapter.length()); x >= matchLength; --x) {
                    if (trimmedSeq.length() >= static_cast<size_t>(x)) {
                        std::string adapterEnd = adapter.substr(adapter.length() - x);
                        std::string toMatch = trimmedSeq.substr(0, x);
                        int mismatches = 0;
                        if (maxMismatchCost > 0)
                            mismatches = calculateEditDistance(adapterEnd, toMatch);
                        else
                            for (size_t i = 0; i < toMatch.length(); ++i)
                                mismatches += calculateBaseDistance(toMatch[i], adapterEnd[i]);
                        bool shouldTrim = (maxMismatchCost > 0) ? (mismatches <= maxMismatchCost) : (mismatches <= maxMismatches);
                        if (shouldTrim) {
                            trimmedSeq = trimmedSeq.substr(x);
                            trimmedQual = trimmedQual.substr(x);
                            usedAdapterIndex = adapterIdx;
                            if (recurse) trimmed = true;
                            break;
                        }
                    }
                }
            }
            if (!trimmed && partial) {
                for (int x = matchLength; x < static_cast<int>(adapter.length()); ++x) {
                    if (trimmedSeq.length() >= static_cast<size_t>(x)) {
                        for (size_t y = 0; y <= adapter.length() - x; ++y) {
                            std::string adapterSub = adapter.substr(y, x);
                            std::string toMatch = trimmedSeq.substr(0, x);
                            int mismatches = 0;
                            if (maxMismatchCost > 0)
                                mismatches = calculateEditDistance(adapterSub, toMatch);
                            else
                                for (size_t i = 0; i < toMatch.length(); ++i)
                                    mismatches += calculateBaseDistance(toMatch[i], adapterSub[i]);
                            bool shouldTrim = (maxMismatchCost > 0) ? (mismatches <= maxMismatchCost) : (mismatches <= maxMismatches);
                            if (shouldTrim) {
                                trimmedSeq = trimmedSeq.substr(x);
                                trimmedQual = trimmedQual.substr(x);
                                usedAdapterIndex = adapterIdx;
                                if (recurse) trimmed = true;
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

    std::string reverseString(const std::string& str) {
        std::string reversed = str;
        std::reverse(reversed.begin(), reversed.end());
        return reversed;
    }
};

// ─── Thread worker result ───
struct ThreadResult {
    std::string out1;
    std::string out2;
    std::string singleOut;
    std::vector<AdapterInfo> stats;
};

static void threadWorker(const std::vector<AdapterInfo>& adapterArray,
                         int matchLength, int maxMismatches, int maxMismatchCost,
                         bool recurse, bool partial, int minReadLength,
                         const std::vector<std::string>& readnames1,
                         const std::vector<std::string>& reads1,
                         const std::vector<std::string>& quals1,
                         const std::vector<std::string>& readnames2,
                         const std::vector<std::string>& reads2,
                         const std::vector<std::string>& quals2,
                         size_t start, size_t end, bool isPaired,
                         ThreadResult& result) {
    // Create a clean copy with zero trim counts so each thread only counts its own NEW trims
    std::vector<AdapterInfo> cleanAdapter;
    cleanAdapter.reserve(adapterArray.size());
    for (size_t a = 0; a < adapterArray.size(); ++a)
        cleanAdapter.emplace_back(adapterArray[a].name, adapterArray[a].sequence);
    AdapterTrimmer trimmer(cleanAdapter, matchLength, maxMismatches, maxMismatchCost, recurse, partial);
    for (size_t i = start; i < end; ++i) {
        if (isPaired) {
            auto [trimmedSeq1, trimmedQual1] = trimmer.trimAdapters(reads1[i], quals1[i]);
            auto [trimmedSeq2, trimmedQual2] = trimmer.trimAdapters(reads2[i], quals2[i]);
            if (trimmedSeq1.size() >= static_cast<size_t>(minReadLength) &&
                trimmedSeq2.size() >= static_cast<size_t>(minReadLength)) {
                result.out1 += readnames1[i] + "\n" + trimmedSeq1 + "\n+\n" + trimmedQual1 + "\n";
                result.out2 += readnames2[i] + "\n" + trimmedSeq2 + "\n+\n" + trimmedQual2 + "\n";
            } else if (trimmedSeq1.size() >= static_cast<size_t>(minReadLength)) {
                result.singleOut += readnames1[i] + "\n" + trimmedSeq1 + "\n+\n" + trimmedQual1 + "\n";
            } else if (trimmedSeq2.size() >= static_cast<size_t>(minReadLength)) {
                result.singleOut += readnames2[i] + "\n" + trimmedSeq2 + "\n+\n" + trimmedQual2 + "\n";
            }
        } else {
            auto [trimmedSeq, trimmedQual] = trimmer.trimAdapters(reads1[i], quals1[i]);
            if (trimmedSeq.size() >= static_cast<size_t>(minReadLength))
                result.out1 += readnames1[i] + "\n" + trimmedSeq + "\n+\n" + trimmedQual + "\n";
        }
    }
    result.stats = trimmer.getTrimStats();
}

// ─── Stream readers ───
static bool isGzipFile(const std::string& path) {
    std::ifstream ifs(path, std::ios_base::binary);
    if (!ifs) return false;
    unsigned char buf[2];
    ifs.read(reinterpret_cast<char*>(buf), 2);
    return (ifs.gcount() == 2 && buf[0] == 0x1f && buf[1] == 0x8b);
}

static size_t readGzipChunk(gzFile& gz, size_t maxRecords,
                            std::vector<std::string>& readnames,
                            std::vector<std::string>& reads,
                            std::vector<std::string>& quals) {
    char buf[4096];
    size_t count = 0;
    auto trimCR = [](std::string& s) {
        if (!s.empty() && s.back() == '\r') s.pop_back();
    };
    while (count < maxRecords) {
        std::string line, sequence, plus, quality;
        if (gzgets(gz, buf, sizeof(buf)) == nullptr) break;
        line = buf;
        if (!line.empty() && line.back() == '\n') line.pop_back();
        trimCR(line);
        if (line.empty() || line[0] != '@') break;
        if (gzgets(gz, buf, sizeof(buf)) == nullptr) break;
        sequence = buf;
        if (!sequence.empty() && sequence.back() == '\n') sequence.pop_back();
        trimCR(sequence);
        if (gzgets(gz, buf, sizeof(buf)) == nullptr) break;
        plus = buf;
        if (!plus.empty() && plus.back() == '\n') plus.pop_back();
        trimCR(plus);
        if (plus.empty() || plus[0] != '+') break;
        if (gzgets(gz, buf, sizeof(buf)) == nullptr) break;
        quality = buf;
        if (!quality.empty() && quality.back() == '\n') quality.pop_back();
        trimCR(quality);
        if (sequence.length() != quality.length()) break;
        readnames.push_back(line);
        reads.push_back(sequence);
        quals.push_back(quality);
        ++count;
    }
    return count;
}

static size_t readPlainChunk(std::ifstream& ifs, size_t maxRecords,
                             std::vector<std::string>& readnames,
                             std::vector<std::string>& reads,
                             std::vector<std::string>& quals) {
    size_t count = 0;
    std::string line;
    auto trimCR = [](std::string& s) {
        if (!s.empty() && s.back() == '\r') s.pop_back();
    };
    while (count < maxRecords) {
        if (!std::getline(ifs, line)) break;
        trimCR(line);
        if (line.empty() || line[0] != '@') break;
        std::string readname = line;
        if (!std::getline(ifs, line)) break;
        trimCR(line);
        std::string sequence = line;
        if (!std::getline(ifs, line)) break;
        trimCR(line);
        if (line.empty() || line[0] != '+') break;
        if (!std::getline(ifs, line)) break;
        trimCR(line);
        std::string quality = line;
        if (sequence.length() != quality.length()) break;
        readnames.push_back(readname);
        reads.push_back(sequence);
        quals.push_back(quality);
        ++count;
    }
    return count;
}

class FastqProcessor {
private:
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
    bool isPaired;
    int numThreads;

public:
    FastqProcessor(const std::string& out1, const std::string& out2,
                   const std::string& single,
                   const std::vector<AdapterInfo>& adapterInfos,
                   int matchLen, int minLen, int maxMismatchCount, int maxCosts,
                   bool rc, bool p, int threads)
        : outputFile1(out1), outputFile2(out2),
          singleOutputFile(single),
          adapterArray(adapterInfos),
          matchLength(matchLen),
          minReadLength(minLen),
          maxMismatches(maxMismatchCount),
          maxMismatchCost(maxCosts),
          recurse(rc),
          partial(p),
          isPaired(!out2.empty()),
          numThreads(threads) {}

    void processChunk(const std::vector<std::string>& readnames1,
                      const std::vector<std::string>& reads1,
                      const std::vector<std::string>& quals1,
                      const std::vector<std::string>& readnames2,
                      const std::vector<std::string>& reads2,
                      const std::vector<std::string>& quals2) {
        size_t total = readnames1.size();
        if (total == 0) return;
        int numThreadsActual = std::min(numThreads, (int)total);
        if (numThreadsActual < 1) numThreadsActual = 1;
        size_t recordsPerThread = (total + numThreadsActual - 1) / numThreadsActual;
        std::vector<ThreadResult> results(numThreadsActual);
        std::vector<std::thread> threads;
        threads.reserve(numThreadsActual);
        for (int t = 0; t < numThreadsActual; ++t) {
            size_t start = t * recordsPerThread;
            size_t end = std::min(start + recordsPerThread, total);
            if (start >= end) break;
            threads.emplace_back(threadWorker,
                std::cref(adapterArray),
                matchLength, maxMismatches, maxMismatchCost,
                recurse, partial, minReadLength,
                std::cref(readnames1), std::cref(reads1), std::cref(quals1),
                std::cref(readnames2), std::cref(reads2), std::cref(quals2),
                start, end, isPaired,
                std::ref(results[t]));
        }
        for (auto& t : threads) t.join();

        std::ofstream out1File(outputFile1, std::ios_base::app);
        std::ofstream out2File;
        std::ofstream singleFile;
        if (isPaired) {
            out2File.open(outputFile2, std::ios_base::app);
            singleFile.open(singleOutputFile, std::ios_base::app);
        }
        for (int t = 0; t < numThreadsActual; ++t) {
            out1File << results[t].out1;
            if (isPaired) {
                out2File << results[t].out2;
                singleFile << results[t].singleOut;
            }
        }

        for (int t = 0; t < numThreadsActual; ++t) {
            for (size_t a = 0; a < adapterArray.size() && a < results[t].stats.size(); ++a)
                adapterArray[a].trimCount += results[t].stats[a].trimCount;
        }
    }

    std::vector<AdapterInfo> getTrimStats() const { return adapterArray; }

    void printTrimStats() const {
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

    void writeTrimStatsToFile(const std::string& filename) const {
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
};

// gzip compression (for output)
namespace gzip_custom {
std::string compress(const std::string& input) {
    if (input.empty()) return "";
    z_stream strm{};
    strm.next_in  = reinterpret_cast<Bytef*>(const_cast<char*>(input.data()));
    strm.avail_in = static_cast<uInt>(input.size());
    if (deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
                     16 + MAX_WBITS, 8, Z_DEFAULT_STRATEGY) != Z_OK)
        throw std::runtime_error("deflateInit2 failed");
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

// Read multi-FASTA adapter file
std::vector<AdapterInfo> readAdaptersFromFASTA(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Could not open adapter file: " + filename);
    std::vector<AdapterInfo> adapters;
    std::string line, currentSequence, currentName;
    bool inSequence = false;
    while (std::getline(file, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!currentSequence.empty() && !currentName.empty()) {
                currentSequence.erase(
                    std::remove_if(currentSequence.begin(), currentSequence.end(), ::isspace),
                    currentSequence.end());
                adapters.emplace_back(currentName, currentSequence);
                currentSequence.clear();
            }
            currentName = line.substr(1);
            size_t spacePos = currentName.find_first_of(" \t");
            if (spacePos != std::string::npos)
                currentName = currentName.substr(0, spacePos);
            inSequence = true;
            continue;
        }
        if (inSequence) currentSequence += line;
    }
    if (!currentSequence.empty() && !currentName.empty()) {
        currentSequence.erase(
            std::remove_if(currentSequence.begin(), currentSequence.end(), ::isspace),
            currentSequence.end());
        adapters.emplace_back(currentName, currentSequence);
    }
    return adapters;
}

// ─── Main ───
int main(int argc, char* argv[]) {
    std::string inputFile1, inputFile2, outputFile1, outputFile2, singleOutputFile, adapterFile;
    int matchLength = 8;
    int minReadLength = 0;
    int maxMismatches = 0;
    int maxMismatchCost = 0;
    bool recurse = false;
    bool partial = false;
    bool gzipout = false;
    int numThreads = 1;
    std::vector<AdapterInfo> adapterArray;

    int opt;
    while ((opt = getopt(argc, argv, "i:I:o:O:s:f:m:l:n:c:t:rpgh")) != -1) {
        switch (opt) {
            case 'i': inputFile1 = optarg; break;
            case 'I': inputFile2 = optarg; break;
            case 'o': outputFile1 = optarg; break;
            case 'O': outputFile2 = optarg; break;
            case 's': singleOutputFile = optarg; break;
            case 'f': adapterFile = optarg; break;
            case 'm': matchLength = std::stoi(optarg); break;
            case 'l': minReadLength = std::stoi(optarg); break;
            case 'n': maxMismatches = std::stoi(optarg); break;
            case 'c': maxMismatchCost = std::stoi(optarg); break;
            case 't': numThreads = std::stoi(optarg); break;
            case 'r': recurse = true; break;
            case 'p': partial = true; break;
            case 'g': gzipout = true; break;
            case 'h':
                std::cout << "Usage: ./trimTSO [args]\n"
                          << "-------required args-------\n"
                          << "-i [input.fastq(.gz)]\n"
                          << "-o [output]\n"
                          << "-f [adapterfile(multiFASTA)]\n\n"
                          << "-------optional args-------\n"
                          << "-I [input_reverse.fastq(.gz)] (paired-end)\n"
                          << "-O [reverse_output]\n"
                          << "-s [single_output]\n"
                          << "-m min_match_length (default: 8)\n"
                          << "-l min_read_length (default: 0)\n"
                          << "-n max_mismatches_count (default: 0)\n"
                          << "-c max_mismatch_cost (default: 0)\n"
                          << "-t num_threads (default: 1)\n"
                          << "-r trim recursively (default: false)\n"
                          << "-p trim partial adapters (default: false)\n"
                          << "-g gzip output (default: false)\n";
                return 1;
            default:
                std::cerr << "Unknown option: " << char(opt) << "\n";
                return 1;
        }
    }

    if (inputFile1.empty() || outputFile1.empty() || adapterFile.empty()) {
        std::cerr << "Error: Missing required arguments (-i, -o, -f)" << std::endl;
        return 1;
    }
    if (!inputFile2.empty() && outputFile2.empty()) {
        std::cerr << "Error: -O is required when -I is specified" << std::endl;
        return 1;
    }
    if (!inputFile2.empty() && singleOutputFile.empty()) {
        std::cerr << "Error: -s is required for paired-end data" << std::endl;
        return 1;
    }

    outputFile1 = outputFile1 + ".fastq";
    if (!outputFile2.empty()) {
        outputFile2 = outputFile2 + ".fastq";
        singleOutputFile = singleOutputFile + ".fastq";
    }

    std::cout << "Configuration:\n"
              << "Input file 1: " << inputFile1 << "\n";
    if (!inputFile2.empty())
        std::cout << "Input file 2: " << inputFile2 << "\n"
                  << "Output file 2: " << outputFile2 << "\n"
                  << "Single output: " << singleOutputFile << "\n";
    std::cout << "Output file 1: " << outputFile1 << "\n"
              << "Adapter file: " << adapterFile << "\n"
              << "Match length: " << matchLength << "\n"
              << "Min read length: " << minReadLength << "\n";
    if (maxMismatchCost > 0)
        std::cout << "Max mismatch cost: " << maxMismatchCost << "\n";
    else
        std::cout << "Max mismatches: " << maxMismatches << "\n";
    std::cout << "Recursive: " << (recurse ? "true" : "false") << "\n"
              << "Partial: " << (partial ? "true" : "false") << "\n"
              << "Gzip output: " << (gzipout ? "true" : "false") << "\n"
              << "Threads: " << numThreads << "\n";

    try { adapterArray = readAdaptersFromFASTA(adapterFile); }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    if (adapterArray.empty()) {
        std::cerr << "Error: No adapter sequences found." << std::endl;
        return 1;
    }
    std::cout << "Loaded " << adapterArray.size() << " adapter sequences:" << std::endl;
    for (const auto& adapter : adapterArray)
        std::cout << adapter.name << ": " << adapter.sequence << std::endl;

    bool input1IsGzip = isGzipFile(inputFile1);
    bool input2IsGzip = inputFile2.empty() ? false : isGzipFile(inputFile2);

    gzFile gz1 = nullptr, gz2 = nullptr;
    std::ifstream ifs1, ifs2;

    if (input1IsGzip) {
        gz1 = gzopen(inputFile1.c_str(), "rb");
        if (!gz1) { std::cerr << "Error: Could not open gzip: " << inputFile1 << std::endl; return 1; }
    } else {
        ifs1.open(inputFile1, std::ios_base::binary);
        if (!ifs1) { std::cerr << "Error: Could not open: " << inputFile1 << std::endl; return 1; }
    }
    if (!inputFile2.empty()) {
        if (input2IsGzip) {
            gz2 = gzopen(inputFile2.c_str(), "rb");
            if (!gz2) { std::cerr << "Error: Could not open gzip: " << inputFile2 << std::endl; return 1; }
        } else {
            ifs2.open(inputFile2, std::ios_base::binary);
            if (!ifs2) { std::cerr << "Error: Could not open: " << inputFile2 << std::endl; return 1; }
        }
    }

    FastqProcessor processor(outputFile1, outputFile2, singleOutputFile,
                             adapterArray, matchLength, minReadLength,
                             maxMismatches, maxMismatchCost, recurse, partial, numThreads);

    size_t recordsPerChunk = static_cast<size_t>(numThreads) * 100000ULL;
    bool isPaired = !inputFile2.empty();

    std::cout << "Processing (records per chunk: " << recordsPerChunk << ")..." << std::endl;
    // Remove existing output files before starting
    auto removeIfExists = [](const std::string& path) {
        if (std::filesystem::exists(path)) {
            std::filesystem::remove(path);
            std::cout << "Removed existing file: " << path << std::endl;
        }
    };
    removeIfExists(outputFile1);
    if (isPaired) {
        removeIfExists(outputFile2);
        removeIfExists(singleOutputFile);
    }
    std::string statsFile = outputFile1.substr(0, outputFile1.find_last_of('.')) + "_trimming_stats.txt";
    removeIfExists(statsFile);
    std::vector<std::string> chunkReadnames1, chunkReads1, chunkQuals1;
    std::vector<std::string> chunkReadnames2, chunkReads2, chunkQuals2;
    size_t totalProcessed = 0;

    while (true) {
        chunkReadnames1.clear(); chunkReads1.clear(); chunkQuals1.clear();
        chunkReadnames2.clear(); chunkReads2.clear(); chunkQuals2.clear();

        size_t n1 = input1IsGzip
            ? readGzipChunk(gz1, recordsPerChunk, chunkReadnames1, chunkReads1, chunkQuals1)
            : readPlainChunk(ifs1, recordsPerChunk, chunkReadnames1, chunkReads1, chunkQuals1);
        if (n1 == 0) break;

        if (isPaired) {
            size_t n2 = input2IsGzip
                ? readGzipChunk(gz2, n1, chunkReadnames2, chunkReads2, chunkQuals2)
                : readPlainChunk(ifs2, n1, chunkReadnames2, chunkReads2, chunkQuals2);
            if (n2 != n1) {
                std::cerr << "Warning: Mismatched count in pair-end chunk (" << n1 << " vs " << n2 << ")" << std::endl;
                size_t minN = std::min(n1, n2);
                chunkReadnames1.resize(minN); chunkReads1.resize(minN); chunkQuals1.resize(minN);
                chunkReadnames2.resize(minN); chunkReads2.resize(minN); chunkQuals2.resize(minN);
            }
        }

        processor.processChunk(chunkReadnames1, chunkReads1, chunkQuals1,
                              chunkReadnames2, chunkReads2, chunkQuals2);
        totalProcessed += chunkReadnames1.size();
        std::cout << "\rProcessed " << totalProcessed << " reads" << std::flush;
    }
    std::cout << "\nCompleted" << std::endl;

    if (input1IsGzip) gzclose(gz1);
    if (!inputFile2.empty() && input2IsGzip) gzclose(gz2);

    processor.printTrimStats();

    processor.writeTrimStatsToFile(statsFile);

    auto compressAndWrite = [](const std::string& f) {
        size_t size = std::filesystem::file_size(f);
        std::vector<char> buffer(size);
        std::ifstream ifs(f, std::ios_base::binary);
        ifs.read(buffer.data(), size); ifs.close();
        std::string compressed = gzip_custom::compress(std::string(buffer.data(), size));
        std::ofstream ofs(f + ".gz", std::ios_base::binary);
        ofs << compressed; ofs.close();
        std::filesystem::remove(f);
    };

    if (gzipout) {
        std::cout << "Compressing output files" << std::endl;
        if (!outputFile1.empty() && std::filesystem::exists(outputFile1)) compressAndWrite(outputFile1);
        if (!outputFile2.empty() && std::filesystem::exists(outputFile2)) compressAndWrite(outputFile2);
        if (!singleOutputFile.empty() && std::filesystem::exists(singleOutputFile)) compressAndWrite(singleOutputFile);
    }

    return 0;
}
