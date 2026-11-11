# trimTSO 内部TSOトリム実装 Kanban

> **目的**: リード内部にあるTSO断片（fwd/RC両配向、1〜3mm変異含む）を厳格シード走査で検出し、プライマー性質に基づく配向別破棄（fwd→5'破棄/RC→3'破棄）で挿入を保持する
> **状態**: 🚀 In Progress (2026-11-11 Hermes kanban にカード化済。下記「Hermes カード対応」参照)
> **関連ファイル**: `/home/takahashi/trimTSO/trimTSO.cpp`, `/home/takahashi/trimTSO/README.md`, `~/.hermes/cache/scratch/tsotest/`
> **計画書**: `plans/2026-11-10-internal-tso-trim.md`

---

## カンバン

| 状態 | ID | タスク | 見積もり | 依存 | 備考 |
|------|----|--------|---------|------|------|
| 📋 **Plan** | **1** | Step 0: poscls.py で内部TSO量のベースライン実測（.32、prodサブセット） | 10分 | なし | 結果を計画書のStep 0に追記 |
| 📋 **Plan** | **2** | trimTSO.cpp: `-x` boolean ＋ `-R <比>`（既定0.7）＋ `-F <bp>`（既定18）フラグ追加（getopt L761・変数・入力チェック・threadWorker経由・AdapterTrimmerメンバー） | 15分 | なし | 既定false=無効、無効で現行バイト同一を先に確認 |
| 📋 **Plan** | **3** | trimTSO.cpp: `internalTrim()` 実装（seedLen=clamp(seedRatio×長さ, seedFloor, 長さ)・**seedFloor未満アダプタは全長完全一致のみ（mmCap=0）**、fwd→5'破棄/RC→3'破棄、mmCap=0.1×長さ） | 35分 | 2 | 計画書のコードそのまま。rcAdapterArrayは再導入しない |
| 📋 **Plan** | **4** | trimTSO.cpp: `trimAdapters()` 冒頭への組み込み（端パスの前） | 10分 | 3 | pass1/pass2・trimCount統計は無改変 |
| 📋 **Plan** | **5** | trimTSO.cpp: `internalTrimCount` 統計＋statsファイル行追加 | 10分 | 4 | 列追加でなく行追加（Total行の後） |
| 📋 **Plan** | **6** | 合成バッテリ検証（internal_test.py: frag互換L=5〜30＋1〜3mm変異版＋**fwd版（TSジャンクション型）**＋**短アダプタ版（全長一致経路）**＋既存回帰） | 40分 | 3,4 | `-x`なし=現行バイト一致が回帰の証明 |
| 📋 **Plan** | **7** | 本番A/B（S5+277-N）: 差分全数監査＋ペア通過劣化チェック＋配向別規則の実データ確認 | 45分 | 6 | 偽陽性1本でも出たらシード比/下限/mmCap設計に差し戻し |
| 📋 **Plan** | **8** | README.md 更新（`-x`/`-R`/`-F`仕様・配向別破棄の理由・全長一致経路の偶然一致曲線・限界・実測値） | 10分 | 7 | コミット・pushはリンの明示指示時 |

## 依存関係図

```
2 (-x/-R/-Fフラグ) ──→ 3 (internalTrim) ──→ 4 (組み込み) ──→ 5 (統計) ──→ 6 (バッテリ)
1 (ベースライン) ─────────────────────────────────────→ 7 (本番A/B) ←── 6
7 ──→ 8 (README)
```

## Hermes カード対応 (2026-11-11 作成・planner 起票)

| 計画ID | Hermes ID | 担当 | 計画Step | 依存(Hermes) |
|--------|-----------|------|----------|--------------|
| 1 | t_6a4210b6 | worker | Step 0 ベースライン実測 (.32 poscls) | — |
| 2 | t_c4af30f0 | worker | Step 1 -x/-R/-F フラグ | — |
| 3 | t_4d72b65d | worker | Step 2 internalTrim() | 2 |
| 4 | t_ea74ac58 | worker | Step 3 trimAdapters() 組込 | 3 |
| 5 | t_126aba24 | worker | Step 4 internalTrimCount 統計 | 4 |
| 6 | t_f1c361b4 | worker | Step 5 合成バッテリ+回帰 | 5 |
| — | t_9f5b4eb1 | reviewer | Step 5 検証ゲート (read-only) | 6 |
| 7 | t_2704dd53 | worker | Step 6 本番A/B (S5+277-N 全数監査) | 1(Hermes), 6 |
| 8 | t_9e92269e | worker | Step 7 README 更新 | reviewerゲート, 本番A/B |

- 全カード body 共通規律: git commit/push/merge/checkout/stash 禁止、右山括弧1文字でのペイロード切り詰め対策（リダイレクト・diff/コード貼付禁止、スクリプトファイル化）、検証はインライン python -c ではなくスクリプト/CLI 実行、deny 時の同一コマンド再試行禁止
- T3/T4 は仕様書のコードブロックに右山括弧が多量にあるため「python スクリプトで仕様 md を読靶向先へ挿入」する手順を body に明記
- 277-N 実パス（.32）: /home/ngsuser/Documents/261029Fastq/temp/S20260277-N_S1_L001/S20260277-N_S1_L001{F,R}_fp.fastq.gz（T8 body に反映済）

## 実装メモ

### 判定済みの設計（2026-11-10 Discord確定）
- 処理ルール: **C案v2（配向別破棄）**。fwd断片=挿入は3'側→断片右端で5'側一括破棄。RC断片=挿入は5'側→断片左端で3'側一括破棄。スプライスはキメラ接合の偽リーディングになるため不採用
- ガードは镜像: fwd=「右端以降の3'フランク≥20bp」、RC=「左端より手前の5'フランク≥20bp」。3'末端のfwd断片（実測7〜16本/サンプル）はリードごと捨てず現行どおり残存
- 複数出現は左から順に切り→再走査（最大4周）→「fwdプライマーの右・RCプライマーの左に挟まれた挿入1本」に収束
- シード長: **固定22不採用**。アダプタごとに `clamp(round(seedRatio×長さ), seedFloor, 長さ)`。**seedRatio/seedFloorは引数（-R/-F）で指定・既定=0.7/18**
- **seedFloor未満アダプタはスキップせず全長完全一致のみで切断**（mmCap=0・変異1mmでも検出不可=誤発火抑止の最小構成）
- mmCap: `max(1, round(0.1×アダプタ長))`（30bp→3mm=tails実測の変異幅）
- 感度限界: 完全一致シードなので≤1mmは確実、2mm以上は無変異シード窓がある場合のみ（シードを下げると誤発火が指数増→受け入れる）
- 同設計（端アンカー連鎖）では解けない（frag実測: int_05〜23全残存）→ 全長走査の新Stageが必須
- 内部Stageは端パスの**前に**置く: `-n 1`連鎖の挿入食い破し（int_30型、123bp喫食）を先に止める効果あり
- `-m`/`-n`とは独立（`-m`/`-n`のチューニングはリン担当）

### 実測根拠（本計画の根拠数字）
- 内部断片残存: frag実測 int_05〜int_23 = out=in−12（3'端ニブルのみ）、frag残=True
- int_30 = 3'パス連鎖がI2(80bp)を食い破って到達 → 190→67bp
- survey（RAW 8サンプル）: TSOfwd22 R1=30%/R2=48%、TSOrc22 R1=13%/R2=5.5%
- tails_3p（4サンプル）: motif出現1〜10回/read、連結周期25〜30bp、変異1〜3mm
- 誤発火確率（シード完全一致・リード250bp・配列2・窓全部）: 30bp→seed21 = ~2.4e-9/リード（~0.0007回/サンプル）、24bp→seed18 = ~0.016回/サンプル
- 全長一致経路（seedFloor未満）の偶然一致: 16bp≈0.035回/サンプル、14bp≈0.55回、12bp≈9回（3e5リード想定）→ `-F`設定はこの曲線で判断
- S5既知値: トリム数119,276（-m 5 -n 1 -l 50 -r -p -g -t 4、2026-11-10）

### テスト資材
- 合成: `uni/frag/`（.32）の genfrag.py/chkfrag.py をローカル複製済み（`~/.hermes/cache/scratch/tsodata/uni/`）
- 既存回帰: `~/.hermes/cache/scratch/tsotest/` の diag_in.fastq / con.fastq / ins.fastq
- 本番入力: S5 = `~/.hermes/cache/scratch/tsodata/S20260281-T_S5_L001{F,R}_fp.fastq.gz`、277-N = .32 `/data/results/S20260277-N_S1_L001/`
