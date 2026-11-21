#!/usr/bin/env python3
"""フィクスチャ長さ検証: manifest.tsv と実ファイルのレコード長が一致するか."""
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
fx = os.path.join(here, "fixtures")
rows = []
with open(os.path.join(fx, "manifest.tsv")) as f:
    next(f)
    for line in f:
        c, rid, l1, l2, note = line.rstrip("\n").split("\t")
        rows.append((c, rid, int(l1), int(l2), note))

def seqs(path):
    out = []
    with open(path) as f:
        lines = f.read().splitlines()
    for i in range(0, len(lines), 4):
        out.append((lines[i], lines[i + 1]))
    return out

bad = 0
by_case = {}
for c in "ABCD":
    r1 = seqs(os.path.join(fx, c + "_R1.fq"))
    r2 = seqs(os.path.join(fx, c + "_R2.fq"))
    exp = [r for r in rows if r[0] == c]
    assert len(r1) == len(r2) == len(exp), (c, len(r1), len(r2), len(exp))
    for (h1, s1), (h2, s2), (cc, rid, l1, l2, note) in zip(r1, r2, exp):
        ok = (h1 == "@" + rid + "/1" and h2 == "@" + rid + "/2"
              and len(s1) == l1 and len(s2) == l2)
        if not ok:
            bad += 1
            print("MISMATCH", c, rid, h1, len(s1), l1, h2, len(s2), l2)
    by_case[c] = (len(r1), sorted(set(len(s) for _, s in r1)),
                  sorted(set(len(s) for _, s in r2)))
for c, v in by_case.items():
    print("case", c, "pairs=", v[0], "R1 lens=", v[1], "R2 lens=", v[2])
print("TSO prefix check:", all(
    s.startswith("AAGCAGTGGTATCAACGCAGAGTACATGGG")
    for c in "ABC" for _, s in seqs(os.path.join(fx, c + "_R1.fq"))))
print("verify:", "FAIL bad=%d" % bad if bad else "OK")
sys.exit(1 if bad else 0)
