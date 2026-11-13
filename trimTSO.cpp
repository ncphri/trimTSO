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
#include <cstdlib>
#include <zlib.h>
#include <unistd.h>
#include <zlib.h>
#include <cstring>
#include <cmath>

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
    double maxMismatchValue;   // -n: >=1=許容ミスマッチ数 / <1=窓長に対する割合(floor(窓長*-n))
    int maxMismatchCost;
    bool recurse;
    bool partial;
    bool internalEnable;
    double seedRatio;
    int seedFloor;
    int maxInwardOffset;   // 末端照合の内側オフセット上限(-K、既定0=従来どおり)
    long long internalTrimCount = 0;
    static constexpr int kKeepFlank = 20;
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

    // 末端側/先頭側 スライスの最小ミスマッチ数
    // (相補鎖側の照合は trimAdapters のRCビューパスが担う。実データで5'RC/3'fwdの
    //  交差候補は誤トリムのみ実測(2026-11-10 S5)のため廃止)
    int matchMismatchCount(const std::string& adapterEnd, const std::string& adapterHead,
                           const std::string& toMatch) {
        int best = 1 << 30;
        for (const std::string* cand : {&adapterEnd, &adapterHead}) {
            if (maxMismatchCost > 0) {
                int cost = calculateEditDistance(*cand, toMatch);
                if (cost < best) best = cost;
                if (best == 0) break;
            } else {
                int mismatches = 0;
                for (size_t i = 0; i < toMatch.length(); ++i)
                    mismatches += calculateBaseDistance(toMatch[i], (*cand)[i]);
                if (mismatches < best) best = mismatches;
                if (best == 0) break;
            }
        }
        return best;
    }

    // 許容ミスマッチ数: -n>=1=数, -n<1=窓長に対する割合(floor(窓長*-n))
    int allowedMismatches(int winlen) const {
        return maxMismatchValue < 1.0
            ? (int)(maxMismatchValue * (double)winlen)
            : (int)maxMismatchValue;
    }

    // 内部TSO走査: 全リード長を走査し、アダプタ(A)とそのRC両方の完全一致シードで
    // 内部断片を検出する。プライマー性質に基づく配向別規則:
    //   fwd断片 = 挿入は断片の3'側 → 断片右端eで5'側を一括破棄(3'フランク保持)
    //   RC断片  = 挿入は断片の5'側 → 断片左端sで3'側を一括破棄(5'フランク保持)
    // シード長=clamp(round(seedRatio*アダプタ長), seedFloor, アダプタ長)。
    // アダプタ長<seedFloorのときは全長一致(-n許容)で検出(mmCap延長なし)。
    // 保持フランク<kKeepFlankの発火は触らない(端パスの担当)。
    // 戻り値: 切断した断片数(0=検出なし)
    int internalTrim(std::string& seq, std::string& qual) {
        if (!internalEnable || seq.length() < 40) return 0;
        int frags = 0;
        while (true) {   // 左から1断片切り→再走査(上限なし)。1切断ごとにリードはseedFloor(≥4)bp以上短くなるので必ず停止
            bool cut = false;
            for (size_t adapterIdx = 0; adapterIdx < adapterArray.size() && !cut; ++adapterIdx) {
                const std::string& A = adapterArray[adapterIdx].sequence;
                const int alen = (int)A.length();
                const bool fullOnly = alen < seedFloor;      // 短アダプタは全長完全一致のみ
                const int seedLen = fullOnly ? alen
                    : std::min(alen, std::max(seedFloor, (int)std::lround(alen * seedRatio)));
                const int mmCap = fullOnly ? 0 : std::max(1, (int)std::lround(alen * 0.1));
                const std::string refs[2] = { A, reverseComplement(A) };   // [0]=fwd, [1]=RC
                for (int o = 0; o < 2 && !cut; ++o) {
                    const std::string& ref = refs[o];
                    for (int i = 0; i + seedLen <= alen && !cut; ++i) {   // fullOnlyなら1窓(=全長)
                        // 案B: シード窓ごとにミスマッチ数を数える走査。許容=-n(>=1=数/<1=窓長割合)。
                        // 全長一致経路(短アダプタ)も -n 許容。許容内の窓はミスマッチ最少優先・
                        // 同数は最左の順に試す。
                        const int allow = allowedMismatches(seedLen);
                        const int nWin = (int)seq.length() - seedLen + 1;
                        std::vector<std::pair<int,int>> cand;   // (p, seedMm)
                        for (int mm0 = 0; mm0 <= allow; ++mm0) {
                            for (int p = 0; p < nWin; ++p) {
                                int mmw = 0;
                                for (int k = 0; k < seedLen; ++k) {
                                    if (seq[p + k] != ref[i + k]) {
                                        ++mmw;
                                        if (allow < mmw) break;
                                    }
                                }
                                if (mmw == mm0) cand.emplace_back(p, mm0);
                            }
                        }
                        for (size_t ci = 0; ci < cand.size() && !cut; ++ci) {
                            const size_t p = (size_t)cand[ci].first;
                            // T8b(an1): the mismatch that trips the stop is never
                            // consumed. Extension budget shrinks by the seed mm so the
                            // total mm (seed + extension, trip counted) stays within
                            // (-n value) + mmCap.
                            const int extCap = std::max(0, mmCap - cand[ci].second);
                            if (o == 0) {
                                // fwd: 右端eまで延長(mmCap、fullOnly=0mm) → [0:e)破棄、[e:]保持
                                int mm = 0; size_t e = p + seedLen; int ai = i + seedLen;
                                while (e < seq.length() && ai < alen) {
                                    if (calculateBaseDistance(seq[e], ref[ai]) == 0) { ++e; ++ai; continue; }
                                    if (++mm > extCap) break;
                                    ++e; ++ai;
                                }
                                if (seq.length() - e >= kKeepFlank) {        // 3'フランク最低長
                                    seq = seq.substr(e); qual = qual.substr(e);
                                    cut = true; ++frags;
                                }
                            } else {
                                // RC: 左端sまで延長(mmCap、fullOnly=0mm) → [s:]破棄、[:s]保持
                                int mm = 0; size_t s = p; int ai = i;
                                while (s > 0 && ai > 0) {
                                    if (calculateBaseDistance(seq[s - 1], ref[ai - 1]) == 0) { --s; --ai; continue; }
                                    if (++mm > extCap) break;
                                    --s; --ai;
                                }
                                if (s >= kKeepFlank) {                       // 5'フランク最低長
                                    seq = seq.substr(0, s); qual = qual.substr(0, s);
                                    cut = true; ++frags;
                                }
                            }
                        }
                    }
                }
            }
            if (!cut) break;
        }
        return frags;
    }

public:
    AdapterTrimmer(const std::vector<AdapterInfo>& adapterInfos, int matchLen = 8, double maxMismatchVal = 0, int maxCosts = 0, bool rc = false, bool p = false,
                   bool internalEn = false, double seedR = 0.7, int seedF = 18, int maxInOff = 0)
        : adapterArray(adapterInfos), matchLength(matchLen), maxMismatchValue(maxMismatchVal), maxMismatchCost(maxCosts), recurse(rc), partial(p),
          internalEnable(internalEn), seedRatio(seedR), seedFloor(seedF), maxInwardOffset(maxInOff) {
    }

    std::pair<std::string, std::string> trimAdapters(const std::string& sequence, const std::string& quality) {
        if (sequence.empty() || quality.empty()) return {sequence, quality};
        std::string seq = sequence;
        std::string qual = quality;
        if (internalTrim(seq, qual) > 0) {           // 先に内部断片を配向別に破棄
            std::lock_guard<std::mutex> lock(statsMutex);
            internalTrimCount++;                      // Step 4のカウンタ
        }
        TrimmedRead forwardTrimmed = forwardTrim(seq, qual);
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

    long long getInternalTrimCount() const {
        std::lock_guard<std::mutex> lock(statsMutex);
        return internalTrimCount;
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
        if (internalEnable) std::cout << "Internal TSO trims: " << internalTrimCount << std::endl;
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
        if (internalEnable) outFile << "Internal_TSO_Trimmed_Fragments\t" << internalTrimCount << "\n";
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
                    // -K: アダプタをリード末端にぴったり(off=0)に加え、1〜K塩基内側に
                    // ずらした位置も照合する。外側に数塩基付いたリード
                    // [X..][adapter][挿入]をアダプタごと一括切断できる。off=0優先=従来優先。
                    for (int off = 0; off <= maxInwardOffset; ++off) {
                        if (trimmedSeq.length() >= static_cast<size_t>(off + x)) {
                            std::string adapterEnd = adapter.substr(adapter.length() - x);
                            // 先頭側: adapterの5'末端x塩基 (RCビューでのR先端断片=3'パス側で照合)
                            std::string adapterHead = adapter.substr(0, x);
                            std::string toMatch = trimmedSeq.substr(off, x);
                            int mismatches = matchMismatchCount(adapterEnd, adapterHead, toMatch);
                            bool shouldTrim = (maxMismatchCost > 0) ? (mismatches <= maxMismatchCost) : (mismatches <= allowedMismatches(x));
                            if (shouldTrim) {
                                trimmedSeq = trimmedSeq.substr(off + x);
                                trimmedQual = trimmedQual.substr(off + x);
                                usedAdapterIndex = adapterIdx;
                                if (recurse) trimmed = true;
                                break;
                            }
                        }
                    }
                }
            }
            if (!trimmed && partial) {
                // -r時は連鎖的に食べ続ける(連結体対応)。ガード: 連鎖継続は「0ミスマッチ」または「x >= 2*matchLength」
                bool again = true;
                do {
                again = false;
                trimmed = false;
                // 最長一致優先(降順)で照合
                for (int x = static_cast<int>(adapter.length()) - 1; x >= matchLength; --x) {
                    for (size_t y = 0; y <= adapter.length() - x; ++y) {
                        std::string adapterSub = adapter.substr(y, x);
                        // -K: 部分一致でも同様に内側オフセットを試す(off=0優先=従来優先)
                        for (int off = 0; off <= maxInwardOffset; ++off) {
                            if (trimmedSeq.length() >= static_cast<size_t>(off + x)) {
                                std::string toMatch = trimmedSeq.substr(off, x);
                                // 部分配列は順配向のみ(候補スライスを両引数に渡す)
                                int mismatches = matchMismatchCount(adapterSub, adapterSub, toMatch);
                                bool shouldTrim = (maxMismatchCost > 0) ? (mismatches <= maxMismatchCost) : (mismatches <= allowedMismatches(x));
                                if (shouldTrim) {
                                    trimmedSeq = trimmedSeq.substr(off + x);
                                    trimmedQual = trimmedQual.substr(off + x);
                                    usedAdapterIndex = adapterIdx;
                                    if (recurse) {
                                        trimmed = true;
                                        again = (mismatches == 0) || (x >= 2 * matchLength);
                                    }
                                    break;
                                }
                            }
                        }
                        if (trimmed) break;
                    }
                    if (trimmed) break;
                }
                } while (again && trimmedSeq.length() >= static_cast<size_t>(matchLength));
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
    long long internalTrimCount = 0;
};

static void threadWorker(const std::vector<AdapterInfo>& adapterArray,
                         int matchLength, double maxMismatchVal, int maxMismatchCost,
                         bool recurse, bool partial, int minReadLength,
                         bool internalEnable, double seedRatio, int seedFloor, int maxInOff,
                         const std::vector<std::string>& readnames1,
                         const std::vector<std::string>& reads1,
                         const std::vector<std::string>& quals1,
                         const std::vector<std::string>& readnames2,
                         const std::vector<std::string>& reads2,
                         const std::vector<std::string>& quals2,
                         const std::vector<std::string>& orphanNames,
                         const std::vector<std::string>& orphanReads,
                         const std::vector<std::string>& orphanQuals,
                         size_t start, size_t end, bool isPaired,
                         ThreadResult& result) {
    // Create a clean copy with zero trim counts so each thread only counts its own NEW trims
    std::vector<AdapterInfo> cleanAdapter;
    cleanAdapter.reserve(adapterArray.size());
    for (size_t a = 0; a < adapterArray.size(); ++a)
        cleanAdapter.emplace_back(adapterArray[a].name, adapterArray[a].sequence);
    AdapterTrimmer trimmer(cleanAdapter, matchLength, maxMismatchVal, maxMismatchCost, recurse, partial,
                           internalEnable, seedRatio, seedFloor, maxInOff);
    const size_t pairedN = readnames1.size();
    for (size_t i = start; i < end; ++i) {
        if (i >= pairedN) {
            // Orphaned mate (its partner record was malformed/skipped):
            // treated exactly like a pair whose mate fell below min length —
            // trimmed and written to the single output.
            size_t o = i - pairedN;
            auto [trimmedSeq, trimmedQual] = trimmer.trimAdapters(orphanReads[o], orphanQuals[o]);
            if (trimmedSeq.size() >= static_cast<size_t>(minReadLength))
                result.singleOut += orphanNames[o] + "\n" + trimmedSeq + "\n+\n" + trimmedQual + "\n";
            continue;
        }
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
    result.internalTrimCount = trimmer.getInternalTrimCount();
}

// ─── Stream readers ───
static bool isGzipFile(const std::string& path) {
    std::ifstream ifs(path, std::ios_base::binary);
    if (!ifs) return false;
    unsigned char buf[2];
    ifs.read(reinterpret_cast<char*>(buf), 2);
    return (ifs.gcount() == 2 && buf[0] == 0x1f && buf[1] == 0x8b);
}

// ─── Streaming line reader (fastp-style): buffered gzread/ifstream with NO line-length limit ───
class FastqLineReader {
private:
    static constexpr size_t BUF_SIZE = 1 << 16; // 64KB decompressed read-ahead
    gzFile gz_ = nullptr;
    std::ifstream ifs_;
    std::string pending_;
    std::string pushback_;
    bool hasPushback_ = false;
    bool eof_ = false;

public:
    bool openGz(const std::string& path) {
        gz_ = gzopen(path.c_str(), "rb");
        return gz_ != nullptr;
    }
    bool openPlain(const std::string& path) {
        ifs_.open(path, std::ios_base::binary);
        return ifs_.is_open();
    }
    void close() {
        if (gz_) { gzclose(gz_); gz_ = nullptr; }
        if (ifs_.is_open()) ifs_.close();
    }
    void pushback(const std::string& line) { pushback_ = line; hasPushback_ = true; }
    // Unbounded getline: assembles lines of any length from buffered reads.
    bool getline(std::string& line) {
        if (hasPushback_) { line.swap(pushback_); hasPushback_ = false; return true; }
        while (true) {
            size_t nl = pending_.find('\n');
            if (nl != std::string::npos) {
                line.assign(pending_, 0, nl);
                pending_.erase(0, nl + 1);
                if (!line.empty() && line.back() == '\r') line.pop_back();
                return true;
            }
            if (eof_) {
                if (!pending_.empty()) {
                    line.swap(pending_);
                    if (!line.empty() && line.back() == '\r') line.pop_back();
                    return true;
                }
                return false;
            }
            char chunk[BUF_SIZE];
            if (gz_) {
                int n = gzread(gz_, chunk, static_cast<unsigned>(BUF_SIZE));
                if (n <= 0) { eof_ = true; continue; }
                pending_.append(chunk, static_cast<size_t>(n));
            } else {
                pending_.resize(pending_.size() + BUF_SIZE);
                ifs_.read(&pending_[pending_.size() - BUF_SIZE], BUF_SIZE);
                std::streamsize got = ifs_.gcount();
                pending_.resize(pending_.size() - BUF_SIZE + static_cast<size_t>(got));
                if (got < static_cast<std::streamsize>(BUF_SIZE)) eof_ = true;
            }
        }
    }
};

// A line qualifies as a sequence line if it is non-empty and purely IUPAC bases.
static bool isDnaLine(const std::string& s) {
    if (s.empty()) return false;
    for (char c : s) {
        switch (toupper(static_cast<unsigned char>(c))) {
            case 'A': case 'C': case 'G': case 'T': case 'U': case 'N':
            case 'R': case 'Y': case 'K': case 'M': case 'S': case 'W':
            case 'B': case 'D': case 'H': case 'V': break;
            default: return false;
        }
    }
    return true;
}

// After a malformed record the read position may be out of sync. Scan forward
// for the next plausible record start: a line beginning with '@' whose NEXT
// line is a pure-DNA sequence line. A bare '@'-leading line is NOT trusted as
// a header because quality strings can also start with '@' (Q31+). The
// sequence line is pushed back so the caller reads it next. `candidate`, if
// given, is tried first as the header candidate.
static bool resyncToNextRecord(FastqLineReader& reader, std::string& header,
                               size_t* dropped = nullptr,
                               const std::string* candidate = nullptr) {
    std::string cand;
    bool haveCand = false;
    size_t drop = 0;
    if (candidate && !candidate->empty() && candidate->front() == '@') {
        cand = *candidate;
        haveCand = true;
    }
    std::string line;
    while (true) {
        if (haveCand) {
            if (!reader.getline(line)) { header = cand; if (dropped) *dropped = drop; return true; } // header at EOF
            if (isDnaLine(line)) { header = cand; reader.pushback(line); if (dropped) *dropped = drop; return true; }
            ++drop; // rejected candidate
            haveCand = false;
        }
        if (!reader.getline(line)) { if (dropped) *dropped = drop; return false; }
        if (!line.empty() && line.front() == '@') { cand = line; haveCand = true; }
        else ++drop;
    }
}

// Reads up to maxRecords well-formed records. Malformed records are skipped
// (with a warning) and reading continues at the next plausible record start.
// `skipped` accumulates the number of malformed records dropped.
static size_t readFastqChunk(FastqLineReader& reader, size_t maxRecords,
                             std::vector<std::string>& readnames,
                             std::vector<std::string>& reads,
                             std::vector<std::string>& quals,
                             size_t& skipped) {
    size_t count = 0;
    std::string header, sequence, plus, quality;
    bool haveHeader = false;
    while (count < maxRecords) {
        if (!haveHeader) {
            if (!reader.getline(header)) break;
            if (header.empty() || header[0] != '@') {
                // Out of sync (stray blank/junk line or drift after a skipped
                // record): rescan for the next valid record start.
                size_t dropped = 0;
                if (!resyncToNextRecord(reader, header, &dropped)) break;
                if (dropped > 0) {
                    std::cerr << "\nWarning: lost sync in FASTQ; discarded " << dropped
                              << " unparsable line(s) before the next record." << std::endl;
                    ++skipped;
                }
            }
            haveHeader = true;
        }
        if (!reader.getline(sequence)) break;
        if (!reader.getline(plus)) break;
        if (plus.empty() || plus[0] != '+') {
            std::cerr << "\nWarning: malformed FASTQ at record \"" << header
                      << "\" (missing '+' line); skipping this record." << std::endl;
            ++skipped;
            // Position is untrusted; resync. The '+' slot line may itself be a
            // header (record that lost both seq and '+' lines), so try it.
            size_t dropped = 0;
            if (!resyncToNextRecord(reader, header, &dropped, &plus)) { haveHeader = false; break; }
            if (dropped > 1) {
                std::cerr << "\nWarning: discarded " << dropped
                          << " further unparsable line(s) after record \"" << header << "\"." << std::endl;
            }
            continue; // header is set; read its seq/plus/qual
        }
        if (!reader.getline(quality)) break;
        if (sequence.length() != quality.length()) {
            std::cerr << "\nWarning: malformed FASTQ at record \"" << header
                      << "\" (seq/qual length mismatch); skipping this record." << std::endl;
            ++skipped;
            haveHeader = false; // 4 lines consumed -> likely aligned at next record
            continue;
        }
        readnames.push_back(header);
        reads.push_back(sequence);
        quals.push_back(quality);
        ++count;
        haveHeader = false;
    }
    return count;
}

// ─── Paired-end name alignment ───
// Compare read names ignoring mate-number suffixes: trailing "/1" or "/2",
// or a trailing " 1:" / " 2:" field (common Illumina conventions). No
// allocation: trims both names in place before comparing.
static bool sameBaseName(const std::string& a, const std::string& b) {
    auto trim = [](const std::string& s) -> size_t {
        size_t len = s.size();
        if (len >= 2 && s[len - 2] == '/' && (s[len - 1] == '1' || s[len - 1] == '2'))
            len -= 2;
        else if (len >= 3 && s[len - 2] == ':' && (s[len - 1] == '1' || s[len - 1] == '2') && s[len - 3] == ' ')
            len -= 3;
        return len;
    };
    size_t la = trim(a), lb = trim(b);
    return la == lb && std::equal(a.begin(), a.begin() + la, b.begin());
}

// Find the entry in names[from..to) whose base name equals `name`; npos if none.
static size_t findBaseName(const std::string& name,
                           const std::vector<std::string>& names, size_t from, size_t to) {
    for (size_t k = from; k < to; ++k)
        if (sameBaseName(name, names[k])) return k;
    return std::string::npos;
}

// ─── Streaming FASTQ writer: plain ofstream or gzwrite stream (fastp-style) ───
class FastqWriter {
private:
    bool isGz_ = false;
    gzFile gz_ = nullptr;
    std::ofstream ofs_;
public:
    bool open(const std::string& path, bool gzipOut) {
        isGz_ = gzipOut;
        if (isGz_) {
            gz_ = gzopen(path.c_str(), "wb6"); // level 6, like fastp default
            return gz_ != nullptr;
        }
        ofs_.open(path, std::ios_base::binary);
        return ofs_.is_open();
    }
    void write(const std::string& data) {
        if (data.empty()) return;
        if (isGz_) {
            unsigned offset = 0;
            while (offset < data.size()) {
                int written = gzwrite(gz_, data.data() + offset,
                                      static_cast<unsigned>(data.size() - offset));
                if (written <= 0) {
                    std::cerr << "Error: gzwrite failed" << std::endl;
                    return;
                }
                offset += static_cast<unsigned>(written);
            }
        } else {
            ofs_ << data;
        }
    }
    void close() {
        if (gz_) { gzclose(gz_); gz_ = nullptr; }
        if (ofs_.is_open()) ofs_.close();
    }
    ~FastqWriter() { close(); }
};

class FastqProcessor {
private:
    std::string outputFile1;
    std::string outputFile2;
    std::string singleOutputFile;
    std::vector<AdapterInfo> adapterArray;
    int matchLength;
    int minReadLength;
    double maxMismatchValue;
    int maxMismatchCost;
    bool recurse;
    bool partial;
    bool internalEnable;
    double seedRatio;
    int seedFloor;
    int maxInwardOffset;
    bool isPaired;
    int numThreads;
    bool gzipOut;
    long long internalTrimCount = 0;
    FastqWriter writer1_, writer2_, writerSingle_;
    bool writersOpen_ = false;

    void ensureWritersOpen() {
        if (writersOpen_) return;
        auto fail = [](const std::string& p) {
            std::cerr << "Error: Could not open output file: " << p << std::endl;
            std::exit(1);
        };
        if (!writer1_.open(outputFile1, gzipOut)) fail(outputFile1);
        if (isPaired) {
            if (!writer2_.open(outputFile2, gzipOut)) fail(outputFile2);
            if (!writerSingle_.open(singleOutputFile, gzipOut)) fail(singleOutputFile);
        }
        writersOpen_ = true;
    }

public:
    FastqProcessor(const std::string& out1, const std::string& out2,
                   const std::string& single,
                   const std::vector<AdapterInfo>& adapterInfos,
                   int matchLen, int minLen, double maxMismatchVal, int maxCosts,
                   bool rc, bool p, int threads, bool gzout = false,
                   bool internalEn = false, double seedR = 0.7, int seedF = 18, int maxInOff = 0)
        : outputFile1(out1), outputFile2(out2),
          singleOutputFile(single),
          adapterArray(adapterInfos),
          matchLength(matchLen),
          minReadLength(minLen),
          maxMismatchValue(maxMismatchVal),
          maxMismatchCost(maxCosts),
          recurse(rc),
          partial(p),
          isPaired(!out2.empty()),
          numThreads(threads),
          gzipOut(gzout),
          internalEnable(internalEn), seedRatio(seedR), seedFloor(seedF), maxInwardOffset(maxInOff) {}

    ~FastqProcessor() { closeWriters(); }

    void closeWriters() {
        if (!writersOpen_) return;
        writer1_.close();
        writer2_.close();
        writerSingle_.close();
        writersOpen_ = false;
    }

    void processChunk(const std::vector<std::string>& readnames1,
                      const std::vector<std::string>& reads1,
                      const std::vector<std::string>& quals1,
                      const std::vector<std::string>& readnames2,
                      const std::vector<std::string>& reads2,
                      const std::vector<std::string>& quals2,
                      const std::vector<std::string>& orphanNames = {},
                      const std::vector<std::string>& orphanReads = {},
                      const std::vector<std::string>& orphanQuals = {}) {
        size_t total = readnames1.size() + orphanNames.size();
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
                matchLength, maxMismatchValue, maxMismatchCost,
                recurse, partial, minReadLength,
                internalEnable, seedRatio, seedFloor, maxInwardOffset,
                std::cref(readnames1), std::cref(reads1), std::cref(quals1),
                std::cref(readnames2), std::cref(reads2), std::cref(quals2),
                std::cref(orphanNames), std::cref(orphanReads), std::cref(orphanQuals),
                start, end, isPaired,
                std::ref(results[t]));
        }
        for (auto& t : threads) t.join();

        ensureWritersOpen();
        for (int t = 0; t < numThreadsActual; ++t) {
            writer1_.write(results[t].out1);
            if (isPaired) {
                writer2_.write(results[t].out2);
                writerSingle_.write(results[t].singleOut);
            }
        }

        for (int t = 0; t < numThreadsActual; ++t) {
            for (size_t a = 0; a < adapterArray.size() && a < results[t].stats.size(); ++a)
                adapterArray[a].trimCount += results[t].stats[a].trimCount;
            internalTrimCount += results[t].internalTrimCount;
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
        if (internalEnable) std::cout << "Internal TSO trims: " << internalTrimCount << std::endl;
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
        if (internalEnable) outFile << "Internal_TSO_Trimmed_Fragments\t" << internalTrimCount << "\n";
        outFile.close();
        std::cout << "Trimming statistics written to: " << filename << std::endl;
    }
};

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
    double maxMismatchValue = 0;
    int maxMismatchCost = 0;
    bool recurse = false;
    bool partial = false;
    bool gzipout = false;
    bool internalEnable = false;
    double seedRatio = 0.7;
    int seedFloor = 18;
    int maxInwardOffset = 0;
    int numThreads = 1;
    std::vector<AdapterInfo> adapterArray;

    int opt;
    while ((opt = getopt(argc, argv, "i:I:o:O:s:f:m:l:n:c:t:xR:F:K:rpgh")) != -1) {
        switch (opt) {
            case 'i': inputFile1 = optarg; break;
            case 'I': inputFile2 = optarg; break;
            case 'o': outputFile1 = optarg; break;
            case 'O': outputFile2 = optarg; break;
            case 's': singleOutputFile = optarg; break;
            case 'f': adapterFile = optarg; break;
            case 'm': matchLength = std::stoi(optarg); break;
            case 'l': minReadLength = std::stoi(optarg); break;
            case 'n': maxMismatchValue = std::stod(optarg); break;
            case 'c': maxMismatchCost = std::stoi(optarg); break;
            case 't': numThreads = std::stoi(optarg); break;
            case 'x': internalEnable = true; break;
            case 'R': seedRatio = std::stod(optarg); break;
            case 'F': seedFloor = std::stoi(optarg); break;
            case 'K': maxInwardOffset = std::stoi(optarg); break;
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
                          << "-n mismatch count (>=1) or mismatch rate per window (<1, e.g. 0.05) (default: 0)\n"
                          << "-c max_mismatch_cost (default: 0)\n"
                          << "-K max_inward_offset_bp_for_end_trim (default: 0)\n"
                          << "-t num_threads (default: 1)\n"
                          << "-r trim recursively (default: false)\n"
                          << "-p trim partial adapters (default: false)\n"
                          << "-g gzip output, streamed directly to .fastq.gz (default: false)\n";
                return 1;
            default:
                std::cerr << "Unknown option: " << char(opt) << "\n";
                return 1;
        }
    }

    if (internalEnable && (seedRatio <= 0 || seedRatio > 1.0 || seedFloor < 4)) {
        std::cerr << "Error: -R must be in (0,1] and -F must be >= 4 when -x is enabled" << std::endl;
        return 1;
    }
    if (maxInwardOffset < 0) {
        std::cerr << "Error: -K must be >= 0" << std::endl;
        return 1;
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

    std::string outputExt = gzipout ? ".fastq.gz" : ".fastq";
    std::string statsFile = outputFile1 + "_trimming_stats.txt";
    outputFile1 = outputFile1 + outputExt;
    if (!outputFile2.empty()) {
        outputFile2 = outputFile2 + outputExt;
        singleOutputFile = singleOutputFile + outputExt;
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
        std::cout << "Max mismatches: " << maxMismatchValue << "\n";
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

    FastqLineReader reader1, reader2;
    if (input1IsGzip) {
        if (!reader1.openGz(inputFile1)) { std::cerr << "Error: Could not open gzip: " << inputFile1 << std::endl; return 1; }
    } else {
        if (!reader1.openPlain(inputFile1)) { std::cerr << "Error: Could not open: " << inputFile1 << std::endl; return 1; }
    }
    if (!inputFile2.empty()) {
        if (input2IsGzip) {
            if (!reader2.openGz(inputFile2)) { std::cerr << "Error: Could not open gzip: " << inputFile2 << std::endl; return 1; }
        } else {
            if (!reader2.openPlain(inputFile2)) { std::cerr << "Error: Could not open: " << inputFile2 << std::endl; return 1; }
        }
    }

    FastqProcessor processor(outputFile1, outputFile2, singleOutputFile,
                             adapterArray, matchLength, minReadLength,
                             maxMismatchValue, maxMismatchCost, recurse, partial, numThreads, gzipout, internalEnable, seedRatio, seedFloor, maxInwardOffset);

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
    removeIfExists(statsFile);
    std::vector<std::string> chunkN1, chunkR1, chunkQ1, chunkN2, chunkR2, chunkQ2;
    // Carry-over buffers: records whose mates are not loaded yet (one side may
    // run ahead of the other after malformed records are skipped).
    std::vector<std::string> pendN1, pendR1, pendQ1, pendN2, pendR2, pendQ2;
    size_t totalProcessed = 0;
    size_t totalSkipped = 0;
    size_t totalOrphans = 0;
    size_t orphanWarnings = 0;
    bool eof1 = false, eof2 = false;
    bool positional = false, nameModeChecked = false;

    while (true) {
        chunkN1.clear(); chunkR1.clear(); chunkQ1.clear();
        chunkN2.clear(); chunkR2.clear(); chunkQ2.clear();

        size_t n1 = readFastqChunk(reader1, recordsPerChunk, chunkN1, chunkR1, chunkQ1, totalSkipped);
        if (n1 < recordsPerChunk) eof1 = true;

        if (!isPaired) {
            if (n1 == 0) break; // 0 = true EOF (bad records are skipped internally)
            processor.processChunk(chunkN1, chunkR1, chunkQ1, chunkN2, chunkR2, chunkQ2);
            totalProcessed += n1;
            std::cout << "\rProcessed " << totalProcessed << " reads" << std::flush;
            continue;
        }

        size_t n2 = readFastqChunk(reader2, recordsPerChunk, chunkN2, chunkR2, chunkQ2, totalSkipped);
        if (n2 < recordsPerChunk) eof2 = true;

        pendN1.insert(pendN1.end(), chunkN1.begin(), chunkN1.end());
        pendR1.insert(pendR1.end(), chunkR1.begin(), chunkR1.end());
        pendQ1.insert(pendQ1.end(), chunkQ1.begin(), chunkQ1.end());
        pendN2.insert(pendN2.end(), chunkN2.begin(), chunkN2.end());
        pendR2.insert(pendR2.end(), chunkR2.begin(), chunkR2.end());
        pendQ2.insert(pendQ2.end(), chunkQ2.begin(), chunkQ2.end());

        // One-time check: do the two files share read names? If not, fall back
        // to positional pairing so unrelated naming conventions still work.
        if (!nameModeChecked && pendN1.size() >= 2 && pendN2.size() >= 2) {
            size_t chk = std::min(std::min(pendN1.size(), pendN2.size()), static_cast<size_t>(1000));
            size_t matches = 0;
            for (size_t k = 0; k < chk; ++k)
                if (sameBaseName(pendN1[k], pendN2[k])) ++matches;
            if (matches * 10 < chk) {
                positional = true;
                std::cerr << "\nWarning: read names do not correspond between -i and -I; "
                             "falling back to positional pairing." << std::endl;
            }
            nameModeChecked = true;
        }

        std::vector<std::string> pairN1, pairR1, pairQ1, pairN2, pairR2, pairQ2;
        std::vector<std::string> orphN, orphR, orphQ;

        auto warnOrphan = [&](const std::string& nm) {
            if (orphanWarnings < 10) {
                ++orphanWarnings;
                std::cerr << "\nWarning: mate of \"" << nm
                          << "\" was lost; writing it to single output." << std::endl;
            }
        };
        auto orphanFrom1 = [&](size_t a, size_t b) {
            for (size_t k = a; k < b; ++k) { warnOrphan(pendN1[k]);
                orphN.push_back(std::move(pendN1[k])); orphR.push_back(std::move(pendR1[k])); orphQ.push_back(std::move(pendQ1[k])); }
        };
        auto orphanFrom2 = [&](size_t a, size_t b) {
            for (size_t k = a; k < b; ++k) { warnOrphan(pendN2[k]);
                orphN.push_back(std::move(pendN2[k])); orphR.push_back(std::move(pendR2[k])); orphQ.push_back(std::move(pendQ2[k])); }
        };

        size_t i = 0, j = 0;
        const size_t s1 = pendN1.size(), s2 = pendN2.size();

        if (positional) {
            size_t m = std::min(s1, s2);
            for (size_t k = 0; k < m; ++k) {
                pairN1.push_back(std::move(pendN1[k])); pairR1.push_back(std::move(pendR1[k])); pairQ1.push_back(std::move(pendQ1[k]));
                pairN2.push_back(std::move(pendN2[k])); pairR2.push_back(std::move(pendR2[k])); pairQ2.push_back(std::move(pendQ2[k]));
            }
            i = m; j = m;
            if (eof2) { orphanFrom1(i, s1); i = s1; }
            if (eof1) { orphanFrom2(j, s2); j = s2; }
        } else {
            // Name-based alignment: both files preserve the same record order,
            // so mismatches come from records skipped as malformed on one side.
            while (i < s1 && j < s2) {
                if (sameBaseName(pendN1[i], pendN2[j])) {
                    pairN1.push_back(std::move(pendN1[i])); pairR1.push_back(std::move(pendR1[i])); pairQ1.push_back(std::move(pendQ1[i]));
                    pairN2.push_back(std::move(pendN2[j])); pairR2.push_back(std::move(pendR2[j])); pairQ2.push_back(std::move(pendQ2[j]));
                    ++i; ++j;
                } else {
                    size_t k = findBaseName(pendN1[i], pendN2, j + 1, s2);
                    if (k != std::string::npos) {
                        orphanFrom2(j, k); j = k; // pendN1[i]'s mate is at k; j..k-1 are orphans
                    } else {
                        size_t l = findBaseName(pendN2[j], pendN1, i + 1, s1);
                        if (l != std::string::npos) {
                            orphanFrom1(i, l); i = l;
                        } else if (eof1 && eof2) {
                            orphanFrom1(i, i + 1); ++i; orphanFrom2(j, j + 1); ++j;
                        } else if (eof2) {
                            orphanFrom1(i, i + 1); ++i; // no more R2 data will arrive
                        } else if (eof1) {
                            orphanFrom2(j, j + 1); ++j; // no more R1 data will arrive
                        } else {
                            break; // mates may still arrive in later chunks; wait
                        }
                    }
                }
            }
            if (i < s1 && eof2) { orphanFrom1(i, s1); i = s1; }
            if (j < s2 && eof1) { orphanFrom2(j, s2); j = s2; }
        }

        auto erasePrefix = [](std::vector<std::string>& v, size_t n) {
            v.erase(v.begin(), v.begin() + n);
        };
        erasePrefix(pendN1, i); erasePrefix(pendR1, i); erasePrefix(pendQ1, i);
        erasePrefix(pendN2, j); erasePrefix(pendR2, j); erasePrefix(pendQ2, j);

        totalOrphans += orphN.size();
        processor.processChunk(pairN1, pairR1, pairQ1, pairN2, pairR2, pairQ2, orphN, orphR, orphQ);
        totalProcessed += pairN1.size() + orphN.size();
        std::cout << "\rProcessed " << totalProcessed << " reads" << std::flush;

        if (eof1 && eof2 && pendN1.empty() && pendN2.empty()) break;
        if (n1 == 0 && n2 == 0 && pendN1.empty() && pendN2.empty()) break; // safety
    }
    std::cout << "\nCompleted" << std::endl;
    if (totalSkipped > 0)
        std::cout << "Skipped " << totalSkipped << " malformed record(s); remaining reads were output normally." << std::endl;
    if (totalOrphans > 0)
        std::cout << "Orphaned mates (partner lost) written to single output: " << totalOrphans << std::endl;

    reader1.close();
    reader2.close();
    processor.closeWriters();

    processor.printTrimStats();

    processor.writeTrimStatsToFile(statsFile);

    return 0;
}
