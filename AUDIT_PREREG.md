# AUDIT PREREG - 045 (written before the audit ran)
Weaknesses of the original run: (a) no standard gene finder baseline ("Prodigal not installed here"), (b) one seed, (c) cold test set dominated by one lineage group (Polyploviricotina, 121 of 200 records), (d) F1 inflated by 89-93% coding bases, so MCC is the honest metric.
Same data (SHA256 re-verified), same record selection, splits, windows and test sets as src/run.py (executed from it).
A1 Add pyrodigal (Prodigal port, meta mode) as baseline on the same cold and random test sets; mark bases inside predicted genes. Report F1/MCC. Pre-specified statement: if pyrodigal MCC > Net MCC + 0.05 on cold, "the small RNN loses to a standard gene caller".
A2 Net seed variance: retrain with torch seeds 0, 1, 2 (same 3 epochs, same train set); report cold and random F1/MCC/AUROC mean and sd. Seed 0 is compared with the originally reported 0.962/0.315 (torch version differs; not bit-comparable).
A3 Cold-set dominance: recompute F1/MCC for ORF caller, pyrodigal and Net (seed 0) on the cold test set EXCLUDING the Polyploviricotina records.
Original gates are not re-fit. Any change in a verdict is written to README/RESULTS.
