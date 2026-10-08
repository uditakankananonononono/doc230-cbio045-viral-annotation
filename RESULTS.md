# CBIO045 results (all numbers from results/results.json, produced by src/run.py, log in results/run.log)
Data: first 60 MB of RefSeq viral.1.genomic.gbff.gz (accession-ordered, not random). 2,110 records met the 1-20 kb + CDS filter (prereg said 3,021 had a CDS; the length filter cut it). 72 lineage groups (taxonomy field), 18 held out. Train 400 records, test 200 cold-lineage, 200 random-split.
Caveat: the cold test set is dominated by Polyploviricotina (121 of 200 records). Per-base labels are 89-93% coding, so F1 is inflated; MCC/AUROC are the honest metrics.

| Split | ORF caller (L=300) F1/MCC | 3-mer LR F1/MCC/AUROC | Conv1d+BiGRU F1/MCC/AUROC |
|---|---|---|---|
| cold lineage | 0.987/0.816 | 0.961/0.000/0.563 | 0.962/0.315/0.807 |
| random | 0.964/0.658 | 0.937/0.000/0.538 | 0.939/0.225/0.795 |

Gates: D2 NOT met (net 0.962 vs ORF 0.987). D3 met (gap, random minus cold, is -0.023; net does not degrade on cold lineages) - but this is not a sign of strength, the net is below the ORF caller on both. 3-mer LR predicts all-coding (MCC 0).
Other directions: D4 per-group F1 (net) 0.82 (Caudoviricetes, n=3) to 0.97; D5 GC bins show no trend; D7 ORF threshold: L=300 best (train F1 0.965; 900 worst 0.902); D8 36% of net errors lie within 30 nt of a CDS edge (computed on the first 300 errors per record); D9 overlapping-gene records F1 0.936 vs 0.950; D10 short genomes (<3 kb) worst, net 0.917. D6 window-size ablation not run (stated in prereg).
Limits: 3 epochs on CPU, 400 training records, no tuning; a simple ORF caller beats the small RNN here. Replication-level benchmark, no novelty claim.

## Audit (second pass) - verdict strengthened: the small RNN loses to a standard gene caller
Prereg AUDIT_PREREG.md committed before running (one variable-name fix to audit.py before any result existed). src/audit.py, results/audit.json, results/audit.log. Data SHA256 re-verified; seed 0 reproduced the original exactly (cold F1 0.9616, MCC 0.315, AUROC 0.807).
- A1 pyrodigal (Prodigal port, meta mode) baseline on the same test sets: cold F1 0.990 / MCC 0.856 (ORF caller 0.987 / 0.816; Net 0.962 / 0.315). Random split: pyrodigal 0.973 / 0.758, ORF 0.964 / 0.658, Net 0.939 / 0.225. Pre-specified statement met: pyrodigal MCC exceeds Net by 0.54 on cold lineages, so the Net loses to a standard gene caller.
- A2 Net over seeds 0,1,2: cold F1 0.961 (sd 0.002), MCC 0.248 (sd 0.051), AUROC 0.791 (sd 0.017); random F1 0.937, MCC 0.165 (sd 0.053), AUROC 0.777. Seed 0 is the best of three, so the original MCC 0.315 was optimistic; the seed-mean is 0.248.
- A3 cold set excluding Polyploviricotina (47% of cold bases; 526,644 bases, other groups): Net F1 0.953 / MCC 0.277, ORF 0.984 / 0.794, pyrodigal 0.987 / 0.844. The ranking is the same without the dominant group.
D2 stays "NOT met" and is now clearer: not a close call. D3 ("gap <= 0.05, met") remains technically true but is not evidence of strength. No claim that RNNs are generally worse: this was 3 epochs, 400 training records, one architecture, CPU.
