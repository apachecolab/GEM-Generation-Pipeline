# Change Log

## Software and Setup

1. Replaced IBM CPLEX with Gurobi 13.0.2, for Apple Silicon compatibility.
2. Consolidated CarveMe, cobra, gurobipy, memote, diamond, and GNU parallel into one conda environment, `model_generation`, installed and verified by `setup_environment.sh`.
3. Changed `config.sh` from resolving paths dynamically on every source to a flat file generated once by `setup_environment.sh`.

---

## Draft Model Generation

1. Added a COBRApy conversion step to produce `.mat` files from CarveMe's SBML models, for Apple Silicon MATLAB compatibility.
2. Wrapped CarveMe and the COBRApy conversion into one script, `genome_to_draftmodels.sh`. Genomes must already be present locally.

---

## Gap-Filling

Split into two scripts, `gapfilling_alt_rxnsets.m` and `gapfilling_selection.m`.

### Alternative Reaction Set Enumeration

1. Enumerates every carbon source in the panel, not just ones the organism is experimentally known to grow on.
2. `maxNumAlt` increased from 10 to 25.
3. Added a configurable cap on each MILP solve (`maxSolveTime`), previously unbounded.
4. Added a mode to rank retry carbon sources by continuous growth rate instead of list order, when that data exists.

### Selection

1. Combos are deduplicated by reaction-set identity, so two combos that resolve to the same reactions aren't tested twice.
2. Added a scoring method parameter (`mcc` by default), extensible to others when continuous growth data is available.
3. The final winner is selected by the active scoring method's score, not by raw accuracy.

---

## Curation

1. Split into two stages: the first stops after `patience` consecutive iterations with no organism improving; the second re-seeds each organism from its own best iteration and stops on the first regression.
2. Added `scoringMethod` and `growthDataMode` parameters (`mcc` by default; `all` expands to every method when continuous growth data is available), matching Selection's parameters.
3. Donor specialist selection is ranked by rank-gap when continuous growth data is available, or by fewest new reactions needed otherwise, resolving ties that were previously left to indexing order.
4. Single-reaction fixes are tried exhaustively before falling back to combination search, and are selected the same way (rank-gap or fewest new false positives).
5. Added a `donorDir` parameter: defaults to this run's own gap-filled models, or can point to a directory of external `.mat` models. External donors use predicted-growth-only eligibility and fewest-new-reactions ranking, since no matched experimental data exists for them.
6. Added `donorReactionMode` for external donors: candidate reactions come from one live FBA optimum (`fba`), or from a precomputed essential-reaction set (`essential`) via a companion sweep script.

---

## Final Model Formatting

1. SBML export runs through COBRApy instead of `writeCbModel`'s own SBML writer, whose `OutputSBML` binary has no Apple Silicon build. A throwaway non-`v7.3` copy of each model is written for the handoff, since `scipy.io.loadmat` can't read `v7.3` files either.
2. Every gene is given `sboTerm SBO:0000243` on SBML export, matching MATLAB's own SBML writer's default (COBRApy's writer does not set this automatically).
3. `comps`/`compNames`/`metComps` are recomputed from each metabolite's own ID suffix rather than carried through from the curated model, where gap-filling and curation's additions and removals leave it stale or misaligned even when its length happens to still match.

---

## Model Verification

1. Turned into an actual script, `memote_verification.sh`; previously only documented as a manual command in a README, never run as part of the pipeline itself.
