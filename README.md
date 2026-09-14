*At*-LSPHERE genome-scale metabolic model generation pipeline
========================

This repository contains a pipeline for generating genome-scale metabolic models of phyllosphere bacteria from the *At*-LSPHERE culture collection, from genome sequences through gap-filling, curation, and quality-verified final models.

## Model generation pipeline

Located in the `scripts/` directory, contains scripts for draft model reconstruction, gap-filling, curation, final model formatting, and MEMOTE quality verification. See [RunGuide.md](RunGuide.md) for step-by-step usage and [SoftwareSetup.md](SoftwareSetup.md) for installation.

## Generated models and reports

Located in `data/models/`, contains draft, gap-filled, curated, and final genome-scale models. The models are provided both in `.mat` format (as outputs of the model generation scripts) and in `.sbml` format. MEMOTE quality reports for each final model are located in `data/models/final/reports/`.

## Carbon source screen results, medium composition, and genomes

Located in `data/experimental/` and `data/genomes/`, contains results of the *in vitro* carbon source screen used to curate the genome-scale models, the minimal medium composition for FBA simulations, and the reference genomes used to generate draft models.

## Change history

See [ChangeLog.md](ChangeLog.md) for a record of methodological changes made in this pipeline.
