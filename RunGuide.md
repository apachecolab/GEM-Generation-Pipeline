# Run Guide

**Installation:** see [SoftwareSetup.md](SoftwareSetup.md)

**Environment setup:** from the Pipeline's root directory:

```bash
bash setup_environment.sh
```

---

## Step 1: Generating Draft Models from Genomes

**Input:** `data/genomes/protein_sequences/*.faa`

**Output:**

| Path | Contents |
|---|---|
| `data/models/draft/carveme/*.xml` | Draft models from CarveMe, not yet gap-filled. |
| `data/genomes/diamond_hits/*.tsv` | DIAMOND hit tables from CarveMe's internal homology search. |
| `data/models/draft/cobrapy/xml/*.xml` | Same models, re-exported through COBRApy. |
| `data/models/draft/cobrapy/mat/*.mat` | Same models as `.mat`. Downstream MATLAB steps load these directly, since MATLAB's own SBML reader (`TranslateSBML`) has broken FBC support on Apple Silicon (see Troubleshooting in SoftwareSetup.md). |

Runs in parallel across available cores.

**Usage:**

```bash
bash scripts/draft_model_generation/genome_to_draftmodels.sh [--cores N] [--solver NAME]
```

**Example:**

```bash
bash scripts/draft_model_generation/genome_to_draftmodels.sh --cores 4 --solver gurobi
```

---

## Step 2: Enumerating Alternative Gap-Filling Reaction Sets

**Input:**
- `data/models/draft/cobrapy/mat/*.mat` (from Step 1)
- `data/genomes/reference_key/AtLSPHERE_RefSeq.mat`
- `data/experimental/medium/min_med_CSourceScreen.mat`
- `data/experimental/carbon_source_screening/binarizerd_CSourceScreen_Jun2024.xlsx`

**Parameters:**

| Parameter | Description |
|---|---|
| `forceO2Uptake` | Sets a mandatory oxygen uptake flux, for aerobic gap-filling. |
| `maxNumAlt` | Maximum number of alternative solutions per carbon source. |
| `maxSolveTime` | Seconds per MILP solve; empty means no limit. |
| `retryMode` | Order of carbon sources tried when gap-filling fails and retries with one forced into the medium: `binary` (true-positive list order) or `continuous` (ranked by growth rate). |
| `growthRateFile` | Required if `retryMode = 'continuous'`. |

**Output:**

| Path | Contents |
|---|---|
| `data/gapfilling/alt_rxnsets/AltRxns_<organism>.mat` | Alternative reaction sets per carbon source (`ActRxnsAll`), plus the carbon source list (`nutrients`). |

**Usage:**

Open MATLAB in `scripts/gapfilling/alt_rxnsets/` and run:

```matlab
gapfilling_alt_rxnsets
```

---

## Step 3: Selecting the Best Gap-Filling Combination

**Input:**
- `data/gapfilling/alt_rxnsets/AltRxns_<organism>.mat` (from Step 2)
- `data/genomes/reference_key/AtLSPHERE_RefSeq.mat`
- `data/experimental/medium/min_med_CSourceScreen.mat`
- `data/experimental/carbon_source_screening/binarizerd_CSourceScreen_Jun2024.xlsx`
- `data/models/draft/cobrapy/mat/*.mat`

**Parameters:**

| Parameter | Description |
|---|---|
| `mandatoryGrowthNutrients` | Carbon sources to prefer solutions for, if within `mandatoryFPRCutoff` of the best overall solution. Only applies to `scoringMethod = 'mcc'`. |
| `mandatoryFPRCutoff` | Upper limit on false positive rate, relative to the best solution, for `mandatoryGrowthNutrients` to be preferred. |
| `maxChooseK` | Maximum number of alternative solutions to combine. |
| `maxBestCombos` | Maximum number of combined solutions to test via FBA. |
| `growthDataMode` | `binary` or `continuous`. `continuous` requires both growth data files; `binary` requires only the binary one. |
| `scoringMethod` | `mcc`, `auc`, or (with `growthDataMode = continuous`) `pearson`, `nrmse`, `wnrmse`, `log-nrmse`, `smape`, or `all`. |
| `expMetric` | Only used when `growthDataMode = continuous`: `growthRate`, `maxOD`, or `endpointOD`. |

**Output:**

| Path | Contents |
|---|---|
| `data/models/gapfilled/<method>/<organism>_gapfilled.mat` | Final gapfilled COBRA model. |
| `data/gapfilling/selection_results/<method>/<organism>_selection.mat` | Scoring and prediction details for the winning combination. |

**Usage:**

Open MATLAB in `scripts/gapfilling/selection/` and run:

```matlab
gapfilling_selection
```

---

## Step 4: Curating False Positive/False Negative Predictions

**Input:**
- `data/models/gapfilled/<method>/<organism>_gapfilled.mat` (from Step 3)
- `data/experimental/medium/min_med_CSourceScreen.mat`
- `data/experimental/carbon_source_screening/binarizerd_CSourceScreen_Jun2024.xlsx`
- `<donorDir>/*.mat`, only if `donorDir` is set
- `<donorDir's parent>/essential/*.mat`, only if `donorReactionMode = 'essential'` (from the optional pre-step below)

**Parameters:**

| Parameter | Description |
|---|---|
| `growthDataMode` | `binary` or `continuous`. `continuous` requires both growth data files; `binary` requires only the binary one. |
| `scoringMethod` | `mcc`, `auc`, or (with `growthDataMode = continuous`) `pearson`, `nrmse`, `wnrmse`, `log-nrmse`, `smape`, or `all`. |
| `expMetric` | Only used when `growthDataMode = continuous`: `growthRate`, `maxOD`, or `endpointOD`. |
| `donorDir` | Empty for internal donors (this run's own gap-filled models). Otherwise a directory of external `.mat` donor models. |
| `donorReactionMode` | `fba` or `essential`. Only used when `donorDir` is set. |
| `maxK` | Maximum number of new reactions to combine when fixing a false negative. |
| `maxNumCombos` | Maximum number of reaction combinations to test per fix. |
| `maxIter` | Maximum number of curation iterations, per stage. |
| `patience` | Consecutive non-improving iterations before the first stage stops. |

**Output:**

| Path | Contents |
|---|---|
| `data/models/curated/<method>/<organism>_curated.mat` | Curated COBRA model. |
| `data/gapfilling/curation_results/<method>/modelAccuracy.mat` | Per-organism confusion matrix and metrics, latest iteration. |
| `data/gapfilling/curation_results/<method>/iterationTracking.mat` | Per-metric score, organism x iteration, both stages concatenated. |
| `data/gapfilling/donor_cache/<donorDir basename>/growth.mat` | Cached donor predicted-growth matrix, only if `donorDir` is set. |

**Optional pre-step, only if `donorReactionMode = 'essential'`:**

Precomputes, for every model in `donorDir`, which reactions are essential for growth on each carbon source. Open MATLAB in `scripts/gapfilling/curation/` and run:

```matlab
donor_essentiality_sweep
```

Writes to `<donorDir's parent>/essential/`. Safe to restart; already-swept donors are skipped.

**Usage:**

Open MATLAB in `scripts/gapfilling/curation/` and run:

```matlab
gapfilling_curation
```

---

## Step 5: Final Model Formatting

**Input:**
- `data/models/curated/<method>/<organism>_curated.mat` (from Step 4)
- `scripts/final_formatting/databases/MNXref4.0.mat`
- `scripts/final_formatting/databases/reactionsRecon.mat`

**Parameters:**

| Parameter | Description |
|---|---|
| `method` | Which `data/models/curated/<method>/` subset to format. |

**Output:**

| Path | Contents |
|---|---|
| `data/models/final/mat/<organism>_final.mat` | Annotated COBRA model. |
| `data/models/final/xml/<organism>.xml` | Annotated SBML model. |

**Usage:**

Open MATLAB in `scripts/final_formatting/` and run:

```matlab
final_model_formatting
```

---

## Step 6: Model Quality Verification

**Input:** `data/models/final/xml/*.xml` (from Step 5)

**Output:**

| Path | Contents |
|---|---|
| `data/models/final/reports/<organism>.html` | MEMOTE quality snapshot report. |

Runs in parallel across available cores.

**Usage:**

```bash
bash scripts/verification/memote_verification.sh [--cores N]
```
