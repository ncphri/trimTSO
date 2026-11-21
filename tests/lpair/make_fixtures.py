#!/usr/bin/env python3
"""合成テスト A-D フィクスチャ生成 (仕様: planner_spec_lpair.md §テスト).

決定論的: 固定 seed。Usage: python3 make_fixtures.py [outdir]
出力先: fixtures/A_R1.fq, A_R2.fq, ..., D_R1.fq, D_R2.fq, adapters.fa, manifest.tsv

注: trimTSO.cpp には TSO 配列定数は存在しない(アダプタは -f FASTA 駆動)。
TSO30 正本配列は /tmp/tso.fa 実測 = 過去ワーク gen_synth.py の TSOF と一致。

ペイロード安全規則により FASTA ヘッダ文字は chr(62) で生成する。
"""
import os
import random
import sys

GT = chr(62)  # FASTA/FQ ヘッダ文字 '>'
TSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG"  # 30bp SMART TSO (正本実測)
SEED = 20261121
COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def safe_bio(n, rng):
    """ランダム生物配列。TSO/TSO-RC の 10bp 以上窓がリード両端で偶発発火
    しないことを保証して返す(端トリム/RC パスの誤発火防止)。"""
    while True:
        s = ''.join(rng.choice("ACGT") for _ in range(n))
        if TSO in s or rc(TSO) in s:
            continue
        bad = False
        for view in (s, rc(s)):
            for x in range(10, 31):
                if view.startswith(TSO[-x:]) or view.startswith(TSO[:x]):
                    bad = True
                    break
            if bad:
                break
        if not bad:
            return s


def fq(name, seq):
    # FASTQ ヘッダは '@'。GT (chr(62)) は FASTA 専用。
    return "@" + name + "\n" + seq + "\n+\n" + "I" * len(seq) + "\n"


def main():
    if len(sys.argv) == 2:
        outdir = sys.argv[1]
    else:
        outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(SEED)

    # アダプタ辞書 (実走 arm と同一: /tmp/tso.fa = TSO 単一エントリ)
    with open(os.path.join(outdir, "adapters.fa"), "w") as f:
        f.write(GT + "TSO\n" + TSO + "\n")

    manifest = []  # (case, id, r1_raw, r2_raw, note)

    def emit(case, rid, r1, r2, note):
        r1p = os.path.join(outdir, case + "_R1.fq")
        r2p = os.path.join(outdir, case + "_R2.fq")
        with open(r1p, "a") as f:
            f.write(fq(rid + "/1", r1))
        with open(r2p, "a") as f:
            f.write(fq(rid + "/2", r2))
        manifest.append((case, rid, len(r1), len(r2), note))

    # 古い生ファイル除去 (再実行可能に)
    for c in "ABCD":
        for m in ("R1", "R2"):
            p = os.path.join(outdir, c + "_" + m + ".fq")
            if os.path.exists(p):
                os.remove(p)

    # A: R1=TSO30+挿入25=55bp, R2=125bp / -l 50
    # 仕様本文の「R1=155bp」は (TSO30+挿入25)=55bp・期待「トリム後 25bp」・
    # B/C の記述 (30+45=75 / 30+125=155) と算術的に矛盾するため誤記と確定。
    # テスト意図 (トリム後 25bp < 50 の passthrough) を満たす 55bp を採用。
    a_ins = safe_bio(25, rng)
    emit("A", "read_A_0001", TSO + a_ins, safe_bio(125, rng),
         "new-spec: t1=25bp passthrough, t2=125bp, ts=0, mate_passthrough=1 (raw R1=55bp=TSO30+25; spec 155bp is a typo)")

    # B: R1=75(TSO30+挿入45で45bp), R2=78(挿入48bp) / 両 mate < 50
    b_ins = safe_bio(45, rng)
    emit("B", "read_B_0001", TSO + b_ins, safe_bio(48, rng),
         "new-spec: t1=0, t2=0, ts=0 (黙って破棄), pair_dropped_both_short=1")

    # C: 回帰 単一コピー R1=155(TSO30+挿入125), R2=125 / -l 50
    c_ins = safe_bio(125, rng)
    emit("C", "read_C_0001", TSO + c_ins, safe_bio(125, rng),
         "new-spec: t1=125bp, t2=125bp, ts=0, passthrough=0")

    # D: e9bdff3 回帰 TSO x2/x3/x5/x10 反復 200 本 (各50). R1 bio=131bp 維持
    d_ids = []
    n = 0
    for k in (2, 3, 5, 10):
        for _ in range(50):
            n += 1
            rid = "read_D_%04d" % n
            emit("D", rid, TSO * k + safe_bio(131, rng), safe_bio(125, rng),
                 "x%d; new-spec: t1=131bp, t2=125bp" % k)
            d_ids.append(rid)

    with open(os.path.join(outdir, "manifest.tsv"), "w") as f:
        f.write("case\tid\tr1_raw_bp\tr2_raw_bp\tnote\n")
        for row in manifest:
            f.write("\t".join(str(x) for x in row) + "\n")

    print("fixtures written to", outdir)
    print("cases:", "A=1 B=1 C=1 D=%d" % n, "| total records:", len(manifest))


if __name__ == "__main__":
    main()
