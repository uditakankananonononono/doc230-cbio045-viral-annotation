# doc230-cbio045-viral-annotation
Audited: pyrodigal (standard caller) cold MCC 0.856 vs this RNN 0.25 (3-seed mean); the small RNN loses clearly (RESULTS.md).
10 preregistered directions for CBIO045 (Viral Genome Annotation with RNNs). See PREREG.md, RESULTS.md, src/run.py, results/. Run: python3 src/run.py (needs the data per data_manifest).
Honest outcome: the small RNN did not beat a simple ORF baseline (D2 gate missed).

## Reproduce / readability audit (2026-10-09)
Data: sources and SHA256 are in data_manifest/. Scripts read fixed paths under /tmp (download the files there first, check the SHA256, then run python3 src/run.py and python3 src/audit.py). Not independently re-run from a clean machine as part of this readability check; the numbers in RESULTS.md come from the original runs on 2 CPU / 2 GB. Input is the first 60,000,001 bytes of RefSeq viral.1.genomic.gbff.gz saved as /tmp/vir/part.gz (a partial file, not the whole RefSeq viral set).
