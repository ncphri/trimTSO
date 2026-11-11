Library amplification using SMART-9N causes repeated insertion of TSO adapter subsequences at the read ends. This tool specifically recognizes the adapter sequence from each read end and removes it, including the repeated adapter subsequences.

This tool is also useful for adapter/primer trimming as it recognizes and trims adapter sequences from the read ends.

Compile

g++ -std=c++17 trimTSO.cpp -lz -pthread -O3 -o trimTSO

Usage

./trimTSO -i InputF.fastq(.gz) [-I InputR..fastq(.gz)] -o output_F [-O output_R -s single] -f adapters.fa [-m 8 -n 1 -l 20 -r -g]

Recommended run (with internal TSO trimming enabled):

    ./trimTSO -i R1.fastq.gz -I R2.fastq.gz -o out_F -O out_R -s single -f adapters.fa -m 8 -n 1 -l 50 -r -p -g -x -t 4


Internal TSO trimming (-x)
--------------------------

Adapter-derived fragments that slipped into the interior of a read during amplification (internal chimeras) cannot be removed by the existing end-anchored trimming, which only matches adapter sequences at the read ends. With `-x`, the tool scans the full length of each read and detects both the adapter itself (fwd) and its reverse complement (RC) internally, discarding the fragment wholesale. Without `-x`, the tool performs only end trimming (no internal scan).

### Options

| Option | Default | Description |
|---|---|---|
| `-x` | off | Enable internal TSO trimming (boolean) |
| `-R ratio` | 0.7 | Seed length = round(ratio × adapter length). Range-checked (0 < r ≤ 1) only when `-x` is active |
| `-F bp` | 18 | Seed-length floor (seedFloor). Range-checked (≥ 4) only when `-x` is active |

On `-m`: the code default is 8. With smaller match lengths (e.g. 5), short slices of the adapter (≤ 7 bp) can fire by chance in the partial-matching chain of the end-trim path and nibble into pure-insert ends: a pure-insert anchor read of 80 bp was measured to be reduced to 67 bp under such a setting. `-m 8` eliminates this. `-m` and `-n` are tuning knobs; adjust them as operation parameters rather than changing the code.

### Detection mechanism

- Seed length = clamp(round(`-R` × adapter length), `-F`, adapter length). With a 30 bp adapter and the defaults this is 21 bp.
- The seed scan tolerates mismatches up to the `-n` value (no hard cap). Candidate windows within tolerance are tried fewest-mismatches first, ties broken by the leftmost position.
- From the seed position the fragment is extended toward its remaining end, and the cut is made where the total mismatches in the cut region (seed + extension) stay within `-n` + mmCap. mmCap = max(1, round(0.1 × adapter length)) — 3 for a 30 bp adapter. Mismatches beyond the stopping criterion are neither consumed nor included in the fragment.
- The scan is repeated after each fragment cut until no further fragment is found (no cap on the number of passes; each cut shortens the read by at least the seed floor, so the scan always terminates). Reads shorter than 40 bp are not scanned internally (fragments shorter than 40 bp are out of scope).

### Short adapters (adapter shorter than seedFloor)

Adapters shorter than `-F` are detected by full-length match tolerating up to `-n` mismatches (no mmCap extension). With `-n` ≥ 1 the accidental-match rate of this path rises further than the 0-mm curve below, which therefore acts as a lower bound.

The theoretical number of accidental full-length 0-mm matches grows exponentially with adapter length (assuming 250 bp reads, 3e5 reads per input):

| Adapter length | Accidental matches per input |
|---|---|
| 16 bp | ≈ 0.035 |
| 14 bp | ≈ 0.55 |
| 12 bp | ≈ 9 |

With `-n 2`, a mismatch-tolerant seed fires ≈ 0.61 times per input — the same danger zone as a 16 bp exact match. At `-n 1` this regime is not reached. If `-F` is set too low, short adapters fall into the full-length path; consult the curve above when choosing it.

### Orientation rules

- fwd fragment detected: the fragment's 5' side is discarded wholesale (3' flank kept)
- RC fragment detected: the fragment's 3' side is discarded wholesale (5' flank kept)

Rationale: adapter fragments act as primers, so the slipped-in insert always sits on the fragment's 5'-boundary side. No splice of the junction is performed — that would create false reads of the chimera junction.

A fire that would leave a kept flank shorter than 20 bp is left untouched (deferred to the existing end paths). On paired data the internal trim acts on the affected side only; the mate is handled by the existing end path (fwd fires keep the 3' flank, RC fires keep the 5' flank), so the pair remains a valid FR pair.

### Sensitivity limits (as implemented)

- Any fragment whose total mismatches in the cut region (seed + extension) stay within `-n` + mmCap is detectable (with `-n 1` and a 30 bp adapter: up to 4 in total; the seed window may consume up to `-n`, the extension the remainder). A 1-mm mutation in the central band of the adapter is detected at `-n 1` (the seed scan is `-n`-tolerant).
- Mutant fragments with 2 or more mismatches are detected only when they contain a seed-length window with no more than `-n` mismatches AND the total in the cut region stays within the cap.
- The full-length path for adapters shorter than seedFloor tolerates up to `-n` mismatches, with no mmCap extension.
- Short reads (< 40 bp) and cases that would keep < 20 bp of flank are not handled by the internal stage; they belong to the existing 5'/3' end paths.

### Statistics

With `-x` active, one extra line is appended right after the Total line of the statistics file (an added line, not a column):

    Internal_TSO_Trimmed_Fragments	(number of reads with internally cut fragments)

Without `-x` the line does not appear.
