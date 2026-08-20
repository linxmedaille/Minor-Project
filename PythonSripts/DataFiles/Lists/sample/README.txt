AllTheBacteria manifests
========================

batch_pipeline.py reads two manifest files from this folder:

  file_list.r0.2.v2.tsv       one row per genome assembly
  atb.bakta.r0.2.status.tsv   one row per Bakta annotation, with its QC status

The full manifests are 431 MB and 284 MB, which is too large to submit. The
copies here are a 100-sample subset covering a single batch
(atb.assembly.r0.2.batch.1), so that the pipeline can still be run.

Running batch_pipeline.py with these files downloads roughly 33 MB + 294 MB of
archives from OSF and processes 100 real genomes.

The results reported in the dissertation were produced with the full manifests,
covering 202,387 genomes.

The full files are published by AllTheBacteria:
  https://allthebacteria.readthedocs.io/en/latest/
