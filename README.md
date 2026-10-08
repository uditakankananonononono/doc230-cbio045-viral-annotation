# doc230-cbio045-viral-annotation
Audited: pyrodigal (standard caller) cold MCC 0.856 vs this RNN 0.25 (3-seed mean); the small RNN loses clearly (RESULTS.md).
10 preregistered directions for CBIO045 (Viral Genome Annotation with RNNs). See PREREG.md, RESULTS.md, src/run.py, results/. Run: python3 src/run.py (needs the data per data_manifest).
Honest outcome: the small RNN did not beat a simple ORF baseline (D2 gate missed).
