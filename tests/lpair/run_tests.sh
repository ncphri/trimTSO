#!/usr/bin/env bash
# tests/lpair 共通実行コマンド (T1 産物。後続カード T3/T4/T5 が使用する)
#
# Usage: bash run_tests.sh [TRIMTSO_BIN] [OUTDIR]
#   TRIMTSO_BIN 既定 /tmp/trimTSO_pre (現行 HEAD e9bdff3 = RED ベースライン)
#   OUTDIR      既定 <this dir>/out
#
# 各ケース A-D を実走 arm と同一フラグで実行し、measure.py で
# t1/t2/ts 本数・長さ・stats 値・新カウンタ有無を OUTDIR/results.txt に集約する。
# 新仕様の期待値は各ケースの stdout ログと manifest.tsv の note 列参照:
#   A: t1=25bp x1 / t2=125bp x1 / ts=0 / mate_passthrough=1
#   B: t1=0 / t2=0 / ts=0 (黙って破棄) / pair_dropped_both_short=1
#   C: t1=125bp x1 / t2=125bp x1 / ts=0 / passthrough=0
#   D: t1=131bp x200 / t2=125bp x200 / ts=0 (e9bdff3 回帰)
#
# 注: ペイロード安全規則によりリダイレクト記号は使わない (tee / rm+touch で代替)。
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="${1:-/tmp/trimTSO_pre}"
OUT="${2:-$HERE/out}"
FX="$HERE/fixtures"

if [ ! -x "$BIN" ]; then
  echo "ERROR: binary not found or not executable: $BIN"
  exit 1
fi
if [ ! -f "$FX/manifest.tsv" ]; then
  echo "fixtures missing; generating via make_fixtures.py"
  python3 "$HERE/make_fixtures.py"
fi
mkdir -p "$OUT"
rm -f "$OUT/results.txt"
touch "$OUT/results.txt"

# 実走 arm (S39 tso) と同一フラグ。-l 50 は仕様 §テスト 指定。
FLAGS=(-m 10 -n 0.05 -K 4 -l 50 -r -p -t 4)

for C in A B C D; do
  echo "=== case $C ===" | tee -a "$OUT/results.txt"
  "$BIN" \
    -i "$FX/${C}_R1.fq" -I "$FX/${C}_R2.fq" \
    -o "$OUT/${C}_t1" -O "$OUT/${C}_t2" -s "$OUT/${C}_ts" \
    -f "$FX/adapters.fa" \
    "${FLAGS[@]}" |& tee "$OUT/${C}_stdout.txt" | tail -n 1
  python3 "$HERE/measure.py" "$OUT" "$C" | tee -a "$OUT/results.txt"
done

echo "results: $OUT/results.txt"
