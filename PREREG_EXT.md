# PREREG_EXT: 045 extension D6 window-size ablation (written and committed with src/ext_d6.py BEFORE any D6 score exists)
Date 2026-10-09. Closes the declared-but-not-run D6 (PREREG.md: "D6 window-size ablation not run"). With it the parent has 10 directions run.
Data, record selection, splits, test sets, B1 ORF baseline and the Conv1d+BiGRU net are exactly those of src/run.py (the script executes run.py's setup code; SHA256 of the data re-verified in data_manifest). Only the window length changes (stride = W/2).
Configs, each trained with torch seeds 0, 1, 2 (3 epochs unless stated, Adam 1e-3, batch 32):
 C900: W=900 (the original setting, retrained with the same code so the comparison uses one pipeline).
 C300: W=300, 3 epochs.
 C300s: W=300, 1 epoch. Control: W=300 has about 3x the windows, so C300 takes about 3x the optimizer steps of C900; C300s matches steps to C900 so a W effect is not confounded with step count.
Metric: cold-lineage MCC (honest metric here, F1 is inflated by 89-93% coding bases), seed-mean; also F1, AUROC, and the random-split numbers.
D6 gate: window size matters if |mean cold MCC(C300) - mean cold MCC(C900)| >= 0.05 AND C300s differs from C900 by >= 0.05 in the same direction. If only C300 differs, the result is reported as confounded with training steps. Seed sd is reported; with 3 seeds and the earlier MCC sd of 0.05, a pass is weak evidence.
Reference lines from the audit, unchanged: pyrodigal cold MCC 0.856, ORF caller 0.816. This batch does not change the finding that the small net loses to a standard gene caller.
Limits: 3 seeds, one train set of 400 records, one held-out lineage assignment dominated by Polyploviricotina (121 of 200 cold records), labels are RefSeq annotations, only two window sizes (300 and 900; W above 900 would drop records shorter than W under run.py's windowing).
