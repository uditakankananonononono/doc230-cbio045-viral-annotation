# CBIO045 results (all numbers from results/results.json, produced by src/run.py, log in results/run.log)
Data: first 60 MB of RefSeq viral.1.genomic.gbff.gz (accession-ordered, not random). 2,110 records met the 1-20 kb + CDS filter (prereg said 3,021 had a CDS; the length filter cut it). 72 lineage groups (taxonomy field), 18 held out. Train 400 records, test 200 cold-lineage, 200 random-split.
Caveat: the cold test set is dominated by Polyploviricotina (121 of 200 records). Per-base labels are 89-93% coding, so F1 is inflated; MCC/AUROC are the honest metrics.

| Split | ORF caller (L=300) F1/MCC | 3-mer LR F1/MCC/AUROC | Conv1d+BiGRU F1/MCC/AUROC |
|---|---|---|---|
| cold lineage | 0.987/0.816 | 0.961/0.000/0.563 | 0.962/0.315/0.807 |
| random | 0.964/0.658 | 0.937/0.000/0.538 | 0.939/0.225/0.795 |

Gates: D2 NOT met (net 0.962 vs ORF 0.987). D3 met (gap, random minus cold, is -0.023; net does not degrade on cold lineages) - but this is not a sign of strength, the net is below the ORF caller on both. 3-mer LR predicts all-coding (MCC 0).
Other directions: D4 per-group F1 (net) 0.82 (Caudoviricetes, n=3) to 0.97; D5 GC bins show no trend; D7 ORF threshold: L=300 best (train F1 0.965; 900 worst 0.902); D8 36% of net errors lie within 30 nt of a CDS edge (computed on the first 300 errors per record); D9 overlapping-gene records F1 0.936 vs 0.950; D10 short genomes (<3 kb) worst, net 0.917. D6 window-size ablation: was not run at first; since run as an extension (see Extension D6 below): gate not met, underpowered and confounded with step count, no window effect established or excluded.
Limits: 3 epochs on CPU, 400 training records, no tuning; a simple ORF caller beats the small RNN here. Replication-level benchmark, no novelty claim.

## Audit (second pass) - verdict strengthened: the small RNN loses to a standard gene caller
Prereg AUDIT_PREREG.md committed before running (one variable-name fix to audit.py before any result existed). src/audit.py, results/audit.json, results/audit.log. Data SHA256 re-verified; seed 0 reproduced the original exactly (cold F1 0.9616, MCC 0.315, AUROC 0.807).
- A1 pyrodigal (Prodigal port, meta mode) baseline on the same test sets: cold F1 0.990 / MCC 0.856 (ORF caller 0.987 / 0.816; Net 0.962 / 0.315). Random split: pyrodigal 0.973 / 0.758, ORF 0.964 / 0.658, Net 0.939 / 0.225. Pre-specified statement met: pyrodigal MCC exceeds Net by 0.54 on cold lineages, so the Net loses to a standard gene caller.
- A2 Net over seeds 0,1,2: cold F1 0.961 (sd 0.002), MCC 0.248 (sd 0.051), AUROC 0.791 (sd 0.017); random F1 0.937, MCC 0.165 (sd 0.053), AUROC 0.777. Seed 0 is the best of three, so the original MCC 0.315 was optimistic; the seed-mean is 0.248.
- A3 cold set excluding Polyploviricotina (47% of cold bases; 526,644 bases, other groups): Net F1 0.953 / MCC 0.277, ORF 0.984 / 0.794, pyrodigal 0.987 / 0.844. The ranking is the same without the dominant group.
D2 stays "NOT met" and is now clearer: not a close call. D3 ("gap <= 0.05, met") remains technically true but is not evidence of strength. No claim that RNNs are generally worse: this was 3 epochs, 400 training records, one architecture, CPU.

## Extension D6: window-size ablation (PREREG_EXT.md and src/ext_d6.py committed together at 71944c5 before any score; results/ext_d6.json, ext_d6.log). Wording corrected after gate review (see "Correction" at the end of this section).
Same data, splits, test sets and net as run.py; only the window changes (stride W/2). 3 torch seeds each; seed 0 of C900 reproduced the original cold MCC 0.315 exactly. A first launch was lost mid-run (workspace) and rerun from scratch; its partial output is kept in results/ext_d6_run1_partial.*.
| Config | Windows | Optimizer steps | Cold MCC mean (sample sd, ddof=1) | Cold AUROC (3 seeds) | Random-split MCC mean (sample sd, ddof=1) |
|---|---|---|---|---|---|
| C900 (W=900, 3 ep) | 4,150 | 390 | 0.248 (0.063) | 0.807, 0.800, 0.767 | 0.165 (0.064) |
| C300 (W=300, 3 ep) | 12,864 | 1,206 | 0.350 (0.011) | 0.809, 0.803, 0.795 | 0.313 (0.016) |
| C300s (W=300, 1 ep) | 12,864 | 402 | 0.225 (0.027) | 0.765, 0.780, 0.730 | 0.149 (0.016) |
(The first version of this table gave population sd, ddof=0: 0.051, 0.009, 0.022 for cold MCC. The sample sd above is the one to use with n=3.)
D6 gate as preregistered (both C300 and the step-matched C300s must differ from C900 by >= 0.05 in the same direction): NOT MET by point estimates. C300 - C900 = +0.103; C300s - C900 = -0.023.
Welch tests on cold MCC (n=3 per arm): C300s - C900 = -0.023, 95% CI about [-0.157, +0.111], p = 0.60. C300 - C900 = +0.103, 95% CI about [-0.047, +0.252], p = 0.10. Neither interval excludes zero, and neither excludes a 0.10 shift in either direction.
What this does and does not show: D6 is UNDERPOWERED and confounded. Window size and training-step count are not separated with 3 seeds; no window effect is established, and none is excluded. The two arms are asymmetric controls: C300s matches C900 on optimizer steps (402 vs 390) but sees about one third of the bases (1 epoch of 12,864 windows of 300 vs 3 epochs of 4,150 windows of 900); C300 matches C900 on bases seen but has 3x the steps (1,206 vs 390). Neither arm isolates W. The clean test that was NOT run: C900 trained about 9 epochs (about 1,170 steps) against C300 at 3 epochs. Within W=300, 402 -> 1,206 steps went with cold MCC 0.225 -> 0.350 (a step or epoch effect at W=300 only; it says nothing about W=900).
Limits: 3 seeds, so wide intervals; the ablation covers only the smaller window (W=300 vs the original 900; run.py cannot do W>900, since records shorter than W are dropped under its windowing); one train set of 400 records; cold test set dominated by one lineage group (121 of 200 records). The earlier statement that all 10 directions are run depends on this partial, underpowered D6. Even the best config (cold MCC 0.35) stays far below pyrodigal (0.856) and the ORF caller (0.816), so the main finding is unchanged.
Correction: the first version of this section (commit 96b8682) said the W=300 gain "tracks 3x training steps", "no window effect is detectable" and that the net is undertrained at 3 epochs. Those readings overstated what 3 seeds can show and are withdrawn; the wording above replaces them.