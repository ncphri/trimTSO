## 新規セッション用 ハンドオフプロンプト

# trimTSO 内部TSOトリム実装 — 継続実装

## コンテキスト

- リポジトリ: `/home/takahashi/trimTSO`（単一ファイルC++ツール。HEAD = `2bdef0c` "Drop cross-orientation adapter candidates"、origin同期済み）
- 現行設計: `AdapterTrimmer::trimAdapters()`（trimTSO.cpp L142-161）が2パスで端アンカー切断:
  - pass1 `forwardTrim()`（L198-266）= リード5'先端から照合（A末端側/A先頭側スライス + A_sub部分一致連鎖）
  - pass2 = `reverseComplement()`でRCビューにして同じ`forwardTrim()` → リード3'先端から照合（TSO-RC側）
  - 交差候補（RC@5' / fwd@3'）は 2026-11-10 の本番実測で誤トリムのみと判明し削除済み（コミット `2bdef0c`）
- 実行オプション: `-m 5 -n 1 -l "$LenThresh" -r -p -g`（PaFuM2パイプラインから呼ばれる。`-m`/`-n`のチューニングはリン担当で本計画では触らない）
- アダプタ辞書は複数配列可（`PaFuM2_tso_dict.fa` = 7配列）。TSO = `AAGCAGTGGTATCAACGCAGAGTACATGGG`（30bp）、TSO-RC = `CCCATGTACTCTGCGTTGATACCACTGCTT`（31bp）。プロトコルごとに長さ・本数が変わる前提で設計する

### 実測済みの背景データ（.32 = ngsuser@192.168.1.32、鍵済み）

1. **内部断片は現行では一切トリムされない**（合成バッテリ実測、`uni/frag/`）:
   - `genfrag.py`が作る `I1(80bp)+RC[:L]+I2(80bp)`（内部にRC断片5〜30bp）を fastp→trimTSO に通すと:
     - `int_05`〜`int_23`: **内部RC断片が残存**（out = in−12bp。12bpは3'端の偶然一致ニブル。frag残=True）
     - `int_30`: 3'パスの`-n 1`連鎖が挿入I2(80bp)を食い破って内部アダプタに到達 → 190→67bp（断片は消えるが挿入123bp喫食）
     - `end_*`/`pfx_*`: 正しく処理
   - ローカル再現ファイル: `/home/takahashi/.hermes/cache/scratch/tsodata/uni/`（f.fq, fp_tso.fastq.fastq, chkfrag.py, genfrag.py）
2. **本番データの内部TSO量**（`uni/survey.txt`、8サンプルのRAW R1/R2）:
   - TSOfwd22（リード内どこかにTSO fwd 22mer完全一致）: R1の30%、**R2の48%**
   - TSOrc22: R1の13%、R2の5.5%
   - 3'末端25merの首位は常に `GTACTCTGCGTTGATACCACTGCTT`（=TSO-RC 3'断片）= 端アンカー対象
3. **3'末端のTSO-RC連結鎖**（`uni/tails_3p_stats.txt`、4サンプル合算）:
   - 1リードあたりmotif出現 1〜10回（1:576, 2:267, 3:320, 4:298, 5:285, 6:207, 7:107, 8:47, 9:8, 10:1）
   - 連結周期は25〜30bp（28bp:538, 25bp:415, 30bp:345…）= TSO-RC連結体チェーン
   - 変異は1〜3mmで現れる（例: `TTACTCTGCGTTGATACCACTGCTT`等）
   - 実例: `…CGCGTACTCTGCGTTGATACCACTGCTTAGCAGTGGTATCAACGCAGATCATG…` = 挿入1 + RC断片 + TSOfwd断片 + 挿入2（内部キメラ。両側に生物学配列）
4. README_adapters.md（.32 `analysis/adapters/`）の結論: 「read 内部キメラ(TSO+アダプタ+polyA の連鎖)はアダプタ切断では解けない。消すなら分割/破棄ルールが別途必要」＋「fastpは内部位置で切断して3'側を捨てるためTSO短アンカーを渡すとペア通過−84.6%」→ 内部対応は trimTSO 側に新設計が要る

## 設計判断（リン確認済み・2026-11-10 Discord）

**同設計（端アンカー+連鎖）のままでは解けない。** 現行は5'先端/3'先端からしか照合しないので、両側に挿入を持つ内部断片は誰にも照合されない（frag実測: int_05〜23全残存）。新Stage「全長走査」を端パスの前に1つ追加する。候補照合機構・統計・single/merged出力は全部流用。

### 処理ルール（配向別・プライマー性質に基づく。リン確定）

プライマーは挿入の5'境界を定義する。内部に見つかった断片が「プライマー」として:
- **fwd(TSO)断片を内部に検出** → 挿入は断片の**3'側** → **断片の右端eで5'側を一括破棄**（3'フランク保持）
- **TSO-RC断片を内部に検出** → 挿入は断片の**5'側**（3'末端のリードスループライマー） → **断片の左端sで3'側を一括破棄**（5'フランク保持）

| 案 | 動作 | 判定 |
|----|------|------|
| A | 現状維持（端のみ） | 不採用（frag実測で残存確認済み） |
| B | スプライス結合（アダプタ部を切り取って両フランク接合） | 不採用（キメラ2分子を繋ぐ偽リーディング→マッピング/アセンブリ汚染） |
| **C** | **配向別破棄: fwd→5'破棄 / RC→3'破棄**。ペア整合性が保てる（R1=5'フランク、R2=3'フランクが残る=正常FRペア）。実装は既存トリムの延長で最小 | **採用** |
| D | 両フランク保持（分割して2出力） | 不採用（YAGNI。出力系の大幅改修） |

**ガード（镜像）**: fwd側は「断片右端以降の3'フランク ≥ 20bp」、RC側は「断片左端より手前の5'フランク ≥ 20bp」。3'末端のfwd断片（実測でサンプル毎7〜16本）はリードごと捨てず**現行どおり残存**（端パスの担当領域）。

**複数出現の収束**: 左から順に1断片切り→再走査（最大4周）。結果として「fwdプライマーの右・RCプライマーの左に挟まれた挿入1本」に収束する。

### シード長（プロトコル非依存。固定22は不採用・リン確定。係数・下限は引数指定）

アダプタごとに個別計算（辞書は複数配列・長さ違いを前提）。**係数と下限は `-R`/`-F` 引数で指定（既定=現在値）**:

```
seedLen = clamp(round(seedRatio × アダプタ長), seedFloor, アダプタ長)   既定: seedRatio=0.7, seedFloor=18
mmCap(境界延長のミスマッチ許容) = max(1, round(0.1 × アダプタ長))
アダプタ長 < seedFloor → スキップせず**全長完全一致のみで検出**（mmCap延長なし=変異検出不可・誤発火抑止）
```

- **下限（seedFloor）は偶然一致の下限線として必須**。シード完全一致の誤発火確率（リード250bp・配列2・窓全部）:
  - 30bp→seed21(既定): ~2.4e-9/リード → **サンプルあたり~0.0007回**
  - 24bp→seed18(下限適用): ~5e-8/リード → **~0.016回/サンプル**
  - 16bpに下げると ~0.6回/サンプル = 挿入誤切断の危険域
- **seedFloor未満の全長一致経路の偶然一致はアダプタ長の指数で増える**（リード250bp・配列2・全長1窓）: 16bp≈0.035回/サンプル、14bp≈0.55回、**12bp≈9回**。→ seedFloorの設定はこの曲線で判断（全長一致に落とすアダプタは実用上≥16bp推奨）
- 感度の限界（既知として明記する）: 完全一致シードなので**≤1mm変異は確実に検出**、2mm以上は「断片内にシード長の無変異窓が存在する場合」のみ。seedFloor未満の全長一致は**0mmのみ**（1mmでも残存）。2mm以上の残存はシードを下げても誤発火が指数で増えるため受け入れる（実測の残存0.1%クラスと同じ枠）

## 全体フロー（新規）

```
fastq入力（-i/-I）
    ↓
[NEW] internalTrim(): 全長走査（アダプタごと: seedLen=clamp(seedRatio×長さ, seedFloor, 長さ)
      ※-R/-F引数指定可、既定0.7/18。AとそのRC両方、完全一致シード+mmCap延長で断片の端を決める。
       seedFloor未満アダプタは全長完全一致のみ）
      fwd断片 → 断片右端から5'側を一括破棄（3'フランク保持）
      RC断片  → 断片左端から3'側を一括破棄（5'フランク保持）
      ※ガード: 保持フランク < 20bp は触らない（端パスの担当）
    ↓
pass1 forwardTrim()（既存・5'先端）
    ↓
pass2 RCビュー forwardTrim()（既存・3'先端）
    ↓
-l 長さフィルタ / single(merged)出力（既存）
    ↓
trimming_stats.txt（既存+内部トリム行）
```

内部Stageを**端パスの前に**置く理由: `-m 5 -n 1`連鎖が挿入を食い破って内部アダプタに到達する前に、厳格一致で先に切れるので int_30 型の挿入123bp喫食が起こらなくなる（frag実測の再発防止）。

## 改修手順

### Step 0 — 内部TSO量のベースライン実測（実装前に現状数字を固める）

**ファイル**: 実行は .32 上。スクリプトは既存 `~/Documents/261029Fastq/analysis/adapters/uni/poscls.py`（fwd/RC別に pos0/内部(深い)/3'末端61-120/3'末端60内 を分類する既存分類器）

**コマンド**:

```bash
# .32上で。prodサブセット（/tmp/prodsub_*/F.fq が無ければ res6 の prep.sh で再作成）
python3 ~/Documents/261029Fastq/analysis/adapters/uni/poscls.py /tmp/prodsub_N1/F.fq /tmp/prodsub_U3/F.fq /tmp/prodsub_NC48/F.fq /tmp/prodsub_S4/F.fq
```

**期待出力**: サンプル別に `fwd 内部(深い)=N(read M)` `RC 内部(深い)=…` の行。この数字が実装後の効果検証の基準になる。**結果を計画のここに追記すること（実測値のまま残す）**。

**Step 0 実測（2026-11-11、.32・kanban t_6a4210b6）**: 4ファイルとも欠けていたため再作成。uni 配下の prep.sh は `/tmp/adp/<tag>/sub_R1.fq` を作る別スクリプト（参照先 `/data/` は .32 に現存せず）だったので、prodsub の元作成元である psurvey.sh と同一条件（`results/<sample>/<sample>F_fp_t_h_sana.fastq.gz` を head -n 1200000 で解凍）で再作成した。read数は prod_survey.txt の元記録と全一致（N1=218251 / U3=33244 / NC48=78680 / S4=34620）。poscls.py 実測（poscls は全ファイル "F.fq" と表示するためサンプル名は n で対応付けた）:

```
N1   (n=218251)  fwd pos0(5′末端)=0(read 0)   内部(深い)=35(read 35)   3′末端61-120bp=30(read 30)   3′末端60bp内=9(read 9)
N1   (n=218251)  RC  pos0(5′末端)=0(read 0)   内部(深い)=32(read 32)   3′末端61-120bp=144(read 143)  3′末端60bp内=52(read 52)
U3   (n=33244)   fwd pos0(5′末端)=0(read 0)   内部(深い)=45(read 35)   3′末端61-120bp=49(read 42)    3′末端60bp内=27(read 30)
U3   (n=33244)   RC  pos0(5′末端)=0(read 0)   内部(深い)=118(read 108) 3′末端61-120bp=52(read 46)    3′末端60bp内=17(read 21)
NC48 (n=78680)   fwd pos0(5′末端)=0(read 0)   内部(深い)=9(read 6)     3′末端61-120bp=17(read 18)    3′末端60bp内=10(read 10)
NC48 (n=78680)   RC  pos0(5′末端)=0(read 0)   内部(深い)=47(read 46)   3′末端61-120bp=42(read 42)    3′末端60bp内=8(read 8)
S4   (n=34620)   fwd pos0(5′末端)=0(read 0)   内部(深い)=92(read 89)   3′末端61-120bp=35(read 35)    3′末端60bp内=12(read 12)
S4   (n=34620)   RC  pos0(5′末端)=0(read 0)   内部(深い)=26(read 26)   3′末端61-120bp=147(read 145)  3′末端60bp内=53(read 53)
```

効果検証基準（内部(深い) のみ抽出）: fwd 内部(深い)= N1:35 / U3:45 / NC48:9 / S4:92、RC 内部(深い)= N1:32 / U3:118 / NC48:47 / S4:26。実装後（`-x` あり）はこの数が 0〜激減していること、3'末端系（60bp内・61-120bp）は不変であること（端パス担当なので触らない）を Step 6 の検証基準にする。

### Step 1 — CLIフラグ `-x`（有効化）＋ `-R`（シード比）＋ `-F`（シード下限）追加

**ファイル**: `/home/takahashi/trimTSO/trimTSO.cpp`
- L761: getopt文字列 `"i:I:o:O:s:f:m:l:n:c:t:rpgh"` → `"i:I:o:O:s:f:m:l:n:c:t:xR:F:rpgh"`
- L762-797 のswitchに追加:
  - `case 'x': internalEnable = true; break;`（値なしboolean・既定false = 現行バイト同一）
  - `case 'R': seedRatio = std::stod(optarg); break;`（シード長のアダプタ長比、既定0.7、範囲(0,1]）
  - `case 'F': seedFloor = std::stoi(optarg); break;`（シード下限bp、既定18、範囲≥4）
- L749付近の変数宣言: `bool internalEnable = false; double seedRatio = 0.7; int seedFloor = 18;`（既定値=現在の計画値）
- 入力チェック（オプション解析直後）: `internalEnable && (seedRatio <= 0 || seedRatio > 1.0 || seedFloor < 4)` なら使用法エラーで終了
- `threadWorker()`（L283付近）のシグネチャに `bool internalEnable, double seedRatio, int seedFloor` を追加し、L302の `AdapterTrimmer` 生成に渡す。threadWorkerの全呼び出し元（main内のスレッド生成箇所）も引数追加
- AdapterTrimmerのメンバー: `bool internalEnable; double seedRatio; int seedFloor;`
- 定数: `static constexpr int kKeepFlank = 20;`（保持フランク最低長）。mmCapは `max(1, round(0.1 × アダプタ長))` を関数内で計算

### Step 2 — internalTrim() 実装（コア）

**ファイル**: `/home/takahashi/trimTSO/trimTSO.cpp`、`AdapterTrimmer` クラス内（L118の`matchMismatchCount`の下あたりに新設）

**コード全体**（そのままコピー可）:

```cpp
    // 内部TSO走査: 全リード長を走査し、アダプタ(A)とそのRC両方の完全一致シードで
    // 内部断片を検出する。プライマー性質に基づく配向別規則:
    //   fwd断片 = 挿入は断片の3'側 → 断片右端eで5'側を一括破棄(3'フランク保持)
    //   RC断片  = 挿入は断片の5'側 → 断片左端sで3'側を一括破棄(5'フランク保持)
    // シード長=clamp(round(seedRatio*アダプタ長), seedFloor, アダプタ長)。
    // アダプタ長<seedFloorのときは全長完全一致のみで検出(mmCap延長なし=誤発火抑止)。
    // 保持フランク<kKeepFlankの発火は触らない(端パスの担当)。
    // 戻り値: 切断した断片数(0=検出なし)
    int internalTrim(std::string& seq, std::string& qual) const {
        if (!internalEnable || seq.length() < 40) return 0;
        int frags = 0;
        for (int iter = 0; iter < 4; ++iter) {   // 左から1断片切り→再走査。鎖は1回で全域除去
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
                        const std::string seed = ref.substr(i, seedLen);
                        size_t p = seq.find(seed);           // 左端の完全一致シードから
                        while (p != std::string::npos && !cut) {
                            if (o == 0) {
                                // fwd: 右端eまで延長(mmCap、fullOnly=0mm) → [0:e)破棄、[e:]保持
                                int mm = 0; size_t e = p + seedLen; int ai = i + seedLen;
                                while (e < seq.length() && ai < alen) {
                                    if (calculateBaseDistance(seq[e], ref[ai]) == 0) { ++e; ++ai; continue; }
                                    if (++mm > mmCap) break;
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
                                    if (++mm > mmCap) break;
                                    --s; --ai;
                                }
                                if (s >= kKeepFlank) {                       // 5'フランク最低長
                                    seq = seq.substr(0, s); qual = qual.substr(0, s);
                                    cut = true; ++frags;
                                }
                            }
                            p = seq.find(seed, p + 1);
                        }
                    }
                }
            }
            if (!cut) break;
        }
        return frags;
    }
```

**実装上の確認ポイント（実装者が必ず自検する）**:
- `reverseComplement` は既存メンバー（L44）をそのまま使う（rcAdapterArrayは2bdef0cで削除済み。再導入しない）
- fullOnly（アダプタ長<seedFloor）はmmCap=0 → 延長walkは最初のミスマッチで即break → 全長完全一致のみ。1mm変異も検出不可=仕様どおり残存
- IUPACアダプタ（`R`/`Y`等）は `calculateBaseDistance` が既存のまま効くが、シードの`find`は完全一致なのでIUPACを含むアダプタは検出が甘くなる=既知の限界（注意点に記載）
- fwd側の延長オーバーウォーク（mmCapまで挿入側に食い込む）は**保持側の先頭≤mmCap塩基**を失う（RC側は破棄側なので無害）。30bpアダプタなら最大3bp
- シード完全一致の偶然一致確率 = 位置数×窓数×4^-seedLen → 30bpアダプタで~2.4e-9/リード ≪ 1 → 許容一致系の誤発火は原理的に起こらない

### Step 3 — trimAdapters() への組み込み（端パスの前に挿入）

**ファイル**: `/home/takahashi/trimTSO/trimTSO.cpp` L142-161

**変更内容**: 冒頭に内部Stageを1行追加するだけ。

```cpp
    std::pair<std::string, std::string> trimAdapters(const std::string& sequence, const std::string& quality) {
        if (sequence.empty() || quality.empty()) return {sequence, quality};
        std::string seq = sequence;
        std::string qual = quality;
        if (internalTrim(seq, qual) > 0) {           // 先に内部断片を配向別に破棄
            std::lock_guard<std::mutex> lock(statsMutex);
            internalTrimCount++;                      // Step 4のカウンタ
        }
        TrimmedRead forwardTrimmed = forwardTrim(seq, qual);
        // …（既存のpass1/pass2/統計/RC戻しは変更なし。入力を seq/qual に置き換える）
    }
```

既存の`forwardTrim(sequence, quality)`呼び出しを`forwardTrim(seq, qual)`に置き換える箇所だけ変える（L144-150）。RCビュー側・`trimCount`統計・戻し処理は**一切変更しない**。

### Step 4 — 統計カウンタ

**ファイル**: `/home/takahashi/trimTSO/trimTSO.cpp`

- `AdapterTrimmer` に `long long internalTrimCount = 0;`（statsMutexで保護）
- `printTrimStats()`（L168付近）/ `writeTrimStatsToFile()`（L182付近）に1行追加:
  - 標準出力: `Internal TSO trims: N`
  - statsファイル: 既存の `Total` 行の**後**に `Internal_TSO_Trimmed_Fragments\t<数>\n` を追記（既存パーサーはTotal行を見るので列追加より行追加が安全）
- threadWorkerの `result.stats` 合算経路（既存のtrimCount合算と同じ場所）に internalTrimCount の合算を追加

### Step 5 — 合成バッテリ検証

**ファイル**: `/home/takahashi/.hermes/cache/scratch/tsotest/`（ローカル）に新設 `internal_test.py`

**テスト1（frag系・genfrag.py互換）**: `I1(80)+RC[:L]+I2(80)` L=5,8,10,12,14,16,18,20,23,30 と、1〜3mm変異版（RC断片の2〜3ヶ所を置換）。
期待: `-x` あり → 全ケースで `RC断片残存=False` かつ **出力=I1(80bp)のみ**（RC断片とI2は3'破棄）。1mm変異版も全検出、2〜3mm版は「無変異シード長窓があれば検出・なければ残存」（限界どおり）。`-x` なし → 現行とバイト一致（回帰確認）。
**fwd版も追加**: `I1(80)+A[:L]+I2(80)`、`I1(80)+A[15:30]+I2(80)`（TSジャンクション型=TSOの3'側断片）→ 期待: **出力=I2(80bp)のみ**（5'側一括破棄）。

**テスト2（短アダプタ＝seedFloor未満の全長一致経路）**: 辞書に12bp/16bpアダプタを加えた状態で同バッテリを走る。
期待: 全長完全一致（0mm）のみ切断・1mm変異は残存・`-F`を12に下げるとそのアダプタはシード経路（mmCap延長あり）に変わることの確認。偶然一致リスクも実測（rand5000で誤発火=0を確認。12bpアダプタなら理論値~0.15回/5000リードなので0が普通、ただし本番サンプル規模では~9回/サンプルになる点に注意）。

**テスト3（既存バッテリ回帰）**: `diag_in.fastq`（fwd_5p/rc_3p/rc_5p/fwd_3p）、con連結体バッテリ、ins5000を `-x` あり/なしで両走。期待: 5'/3'端ケースは `-x` の有無で**バイト一致**（内部Stageは端アンカーに干渉しないことの証明）。rc_5p/fwd_3pは引き続き残存（正仕様）。

**コマンド**:

```bash
g++ -std=c++17 trimTSO.cpp -lz -pthread -O3 -o trimTSO
cd /home/takahashi/.hermes/cache/scratch/tsotest
./internal_test.py && /home/takahashi/trimTSO/trimTSO -i internal_in.fastq -o internal_out -f adapters.fa -m 5 -n 1 -l 1 -r -p -x
/home/takahashi/trimTSO/trimTSO -i diag_in.fastq -o diag_x -f adapters.fa -m 5 -n 1 -l 1 -r -p -x
/home/takahashi/trimTSO/trimTSO -i internal_in.fastq -o internal_short -f adapters_short.fa -m 5 -n 1 -l 1 -r -p -x -F 18
```

### Step 5 追記 — 本番 `-m` 運用 5→8 へ（リン決定 2026-11-11）

本番運用の `-m` を 5 から 8 に変更する。`-m`/`-n` はチューニング領域であり、コード既定値は変えない（運用フラグ指定として扱う）。以後、検証・本番A/B・README・統計で使う実行例の基準フラグは:

    -m 8 -n 1 -l 50 -r -p -g -x

**根拠（実測+機構分析）**: `-m 5` では 5bp スライスが 1mm 許容で偶然発火し、純挿入端を喫食する。I1 単独（80bp）が 67bp まで削られた（実測。証跡: `tsotest/probe_nibble.py`）。機構はアダプタと部分一致する短スライスが mm 許容連鎖で延長・発火する経路:
- A[19:26] GAGTACA が 1mm で発火（7bp）
- A[9:14] TATCA が 1mm で発火（5bp）
- RC ビューでは同箇所 GGTATC が 1mm で発火（6bp）

`-m 5` ではこれらの 7bp 以下スライスが発火閾値に届くが、`-m 8` では届かず挿入端の喫食が解消される。これに伴い Step 5 のテスト FLAGS・Step 6 の本番引数・Step 7 の README 実行例・T7 の合否基準（旧 67bp 想定は無効、-m 8 再実測が基準）はすべて `-m 8` 運用に統一する（T7/T8/T9 各カードの 2026-11-11 プランナーコメント参照）。

### Step 6 — 本番A/B検証（S5 + 277-N）

**ファイル**: ローカル `/home/takahashi/.hermes/cache/scratch/tsodata/`（S5入力済み）。277-Nは .32 から `S20260277-N_S1_L001_R1_001.fastq.gz`/`_R2_` を取得

**手順**:
1. S5を `-x` あり/なしで両走（`-m 5 -n 1 -l 50 -r -p -g -t 4`、既定の `-R 0.7 -F 18`）
2. 差分リードを全数リストアップし、**全ての内部発火について喫食領域がアダプタ様（AまたはRCに≤mmCap）か目視+スクリプト監査**（本計画の直前実測と同じ手法: full vs 制限版diff→喫食部のアダプタ最短距離）
3. 期待: ①差分リードの喫食部=アダプタ断片のみ（挿入喫食≈0。mmCapまで）、②`-l 50`通過数の落ち込みは内部キメラ分だけ、③トリム数119,276（S5既知値）からの増分=内部トリム数、④fwd発火の出力=3'フランク/RC発火の出力=5'フランク（配向別規則の実データ確認）
4. 277-Nで同様に確認（ペア通過数の劣化が−84.6%級になっていないこと=fastp方式の罠を踏んでいない証明）

### Step 7 — README更新とコミット

**ファイル**: `/home/takahashi/trimTSO/README.md`

- `-x`（boolean・既定無効）、`-R <比>`（既定0.7）、`-F <bp>`（既定18）の説明。seedFloor未満アダプタ=全長完全一致のみ（mmCap延長なし）と偶然一致の曲線（16bp≈0.035/14bp≈0.55/12bp≈9回・サンプルあたり）
- 配向別規則（fwd→5'破棄/RC→3'破棄）とその理由（プライマー=挿入の5'境界。スプライスしない理由=キメラ接合の偽リーディング）
- 感度の限界（≤1mm確実/2mm以上は条件付き、seedFloor未満は0mmのみ）を明記
- 検証結果（Step 5/6の実測値）を追記
- コミット+push（リンの明示指示時に実行）

## 変更ファイル一覧

| ファイル | 変更内容 |
|---------|---------|
| `/home/takahashi/trimTSO/trimTSO.cpp` | `-x`boolean・`-R <比>`・`-F <bp>`フラグ（getopt L761・変数・入力チェック・threadWorker引数）、`internalTrim()`新設（配向別破棄+アダプタ長比例シード+seedFloor未満全長一致経路）、`trimAdapters()`冒頭への組み込み（L142-161）、`internalTrimCount`統計、`printTrimStats`/`writeTrimStatsToFile`の1行追加 |
| `/home/takahashi/trimTSO/README.md` | `-x`/`-R`/`-F`仕様（比例シード+下限・全長一致経路の偶然一致曲線）と配向別規則・限界・検証結果追記 |
| `/home/takahashi/.hermes/cache/scratch/tsotest/internal_test.py` | **新規**: 合成バッテリ（frag互換+fwd版+変異版+短アダプタ版+既存回帰） |

## 注意点

1. **`-x` の既定は無効**。無効なら現行とバイト同一であることをStep 5の回帰テストで必ず証明する（本番パイプラインへの影響ゼロを先に固定してから、パイプライン側の呼び出し変更は別判断にする）
2. **`-m`/`-n` には触れない**: 内部Stageは`-m 5 -n 1`と独立（完全一致シード+アダプタ長比例mmCap）。ただし端パスの連鎖が挿入を食い破る問題（int_30型、123bp喫食）は内部Stageが先に切れることで緩和されるが、`-m 8`等のチューニング判断は引き続きリン担当
3. **統計ファイル形式**: 列追加ではなく行追加（`Internal_TSO_Trimmed_Fragments`はTotal行の後）。PaFuM2側のstatsパーサーがTotal行参照であることを確認（過去ログ解析で確認済みの形式）
4. **ペア整合性**: R1が内部トリムされたペアはR2側は何もしない（R2の端パスは既存どおり）。R1=5'フランク/R2=3'フランクが残るため正常FRペアとして扱える。片方が-l以下になった場合のsingle(merged)振り分けは既存機構がそのまま効く
5. **残存リミット**: 完全一致シードの限界として2mm以上の変異断片は条件付き検出（断片内にシード長の無変異窓がある場合のみ）。seedFloorを下げると誤発火が指数で増える（16bp≈0.6回/サンプル）ので注意
6. **seedFloor未満アダプタの全長一致経路は偶然一致とトレード**: 理論値は16bp≈0.035回/サンプル、14bp≈0.55回、12bp≈9回（リード250bp・3e5リード想定）。`-F`を高くしすぎると短いアダプタがこの経路に入る→曲線を見て設定。この経路はmmCap=0（0mmのみ）なので変異断片は残存=既知の限界としてREADMEに明記
7. **fwd発火で3'末端に届くケースは触らない**: 3'フランク<20bpのfwd断片（3'末端のfwd断片=実測7〜16本/サンプル）はガードがスキップする。リードごと捨てない（fastp方式の−84.6%罠回避）
8. **検証規律**: 本番出力の差分監査は全数（サンプリングしない）。偽陽性が1本でも出たらシード比/下限/mmCapの再設計に戻る

## 実装メモ

### 内部断片の実データ構造（tails_3p_stats.txt）
- 挿入 + [RC断片×n] の3'鎖 = 端アンカー連鎖で既に処理（4230fb7）
- 挿入1 + 鎖 + 挿入2 = 内部キメラ → 今回の対象（RC→3'破棄で挿入1を救済、fwd→5'破棄で挿入2を救済）
- R2は5'側にTSOfwd鎖・3'側にTSO-RC（survey R2 TSOfwd22=48%の正体）→ 端パスで処理済み。R2内部の残り=鎖+挿入の挟み込み=同じStageで処理
- 実例の読み方: `挿入1 + RC断片 + TSOfwd断片 + 挿入2` → 左からfwd断片が先に見つかる→fwd発火=5'破棄→挿入2が残る→再走査→RC断片はすでに消滅→収束。**結果的に一番右のfwdプライマーの右側が残る**

### シードの選び方（プロトコル非依存化の根拠）
- `-R`/`-F` で引数指定（既定=現在値0.7/18）。比率0.7: 30bp→21bp。survey.pyの測定単位（TSOfwd22/TSOrc22）と同程度
- 下限18bp: 偶然一致確率4^-18≈1.45e-11 → サンプルあたり~0.016回（24bpアダプタの場合）。16bpだと~0.6回/サンプル
- seedFloor未満のアダプタは**全長完全一致のみ**（スキップしない・リン確定）。偶然一致はアダプタ長の指数で増える: 16bp≈0.035回/サンプル、14bp≈0.55回、12bp≈9回 → seedFloorの設定はこの曲線で判断
- mmCap=0.1×長さ: tails実測の変異1〜3mm（30bp）に一致。長いアダプタでは比例的に余裕が増えるが、偶然一致の延長オーバーは≤mmCap塩基で収まる
