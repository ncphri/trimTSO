#!/usr/bin/env python3
"""trimTSO 出力 1 ケースの実測: t1/t2/ts の本数・長さ分布 + stats ファイル内容
+ 新カウンタ (mate_passthrough / pair_dropped_both_short) の有無.

Usage: python3 measure.py OUTDIR CASE
出力: 標準出力に 1 ケース分のサマリ (run_tests.sh が results.txt に集約)
"""
import os
import sys
from collections import Counter


def read_fastq(path):
    """素の fastq を全レコード読む。ファイル無しなら None、0 本なら []。"""
    if not os.path.exists(path):
        return None
    seqs = []
    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().strip()
            f.readline()  # '+'
            f.readline()  # qual
            seqs.append(s)
    return seqs


def summarize(recs):
    if recs is None:
        return "MISSING"
    c = Counter(len(s) for s in recs)
    lens = ", ".join("%dbp x%d" % (L, n) for L, n in sorted(c.items()))
    return "n=%d lens=[%s]" % (len(recs), lens)


def main():
    out, case = sys.argv[1], sys.argv[2]
    t1 = read_fastq(os.path.join(out, case + "_t1.fastq"))
    t2 = read_fastq(os.path.join(out, case + "_t2.fastq"))
    ts = read_fastq(os.path.join(out, case + "_ts.fastq"))
    print("[%s] t1 %s | t2 %s | ts %s" % (case, summarize(t1), summarize(t2), summarize(ts)))

    stats_path = os.path.join(out, case + "_t1_trimming_stats.txt")
    blob = ""
    if os.path.exists(stats_path):
        with open(stats_path) as f:
            blob = f.read()
        for line in blob.splitlines():
            print("[%s] stats: %s" % (case, line))
    else:
        print("[%s] stats: MISSING (%s)" % (case, stats_path))
    for key in ("mate_passthrough", "pair_dropped_both_short"):
        print("[%s] counter %s: %s" % (case, key,
              "PRESENT" if key in blob else "ABSENT"))


if __name__ == "__main__":
    main()
