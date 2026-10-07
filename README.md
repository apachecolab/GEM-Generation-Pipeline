Genome-scale metabolic model generation pipeline for environmental bacteria
===========================================================================

---

This repository contains a pipeline for generating genome-scale metabolic models (GEMs) for environmental bacterial isolates, by integrating genome sequences and experimental carbon source utilization data. The pipeline progresses through the generation of draft metabolic reconstructions,<sup>1</sup> gap-filling,<sup>2</sup> curation,<sup>3</sup> and quality verification,<sup>4</sup> to produce curated GEMs that reflect the observed growth characteristics of the organisms on individual carbon sources.

The use of the pipeline is exemplified here through the generation of curated models for 14 strains from the *At*-LSPHERE phyllosphere culture collection,<sup>3,5</sup> including annotated genomes, draft reconstructions, and curated GEMs.

---

## Repository structure

1. **Software setup**\
   Refer to [SoftwareSetup.md](SoftwareSetup.md) for software installation and system setup.

   &nbsp;

2. **Run guide**\
   Refer to [RunGuide.md](RunGuide.md) for detailed instructions and a step-by-step guide to translating genome sequences into curated genome-scale metabolic models.

   &nbsp;

3. **Pipeline inputs**

   a. **Genomes**\
   Located in [`data/genomes/`](data/genomes/). Annotated protein sequences (`.faa`) are in [`protein_sequences/`](data/genomes/protein_sequences/), and [`reference_key/`](data/genomes/reference_key/) maps each RefSeq assembly to its strain name.

   &nbsp;

   b. **Carbon source screen results**\
   Located in [`data/experimental/carbon_source_screening/`](data/experimental/carbon_source_screening/). Binarized growth and no-growth calls for each strain on each carbon source, from the *in vitro* carbon source screen.

   &nbsp;

   c. **Medium composition**\
   Located in [`data/experimental/medium/`](data/experimental/medium/). The minimal medium, vitamins, and carbon sources used for flux balance analysis (FBA) simulations.

   &nbsp;

4. **Capabilities**

   a. **Configurable donor pool for model curation**\
   Curation fixes false negative predictions by adding reactions from donor models. Both parameters are set in [`gapfilling_curation.m`](scripts/gapfilling/curation/gapfilling_curation.m).

   - `donorDir` sets the directory of donor models. By default, these are the gap-filled models of the strains in the same collection.
   - `donorReactionMode` sets which donor reactions are tested: the reactions active in one FBA solution of the donor (`fba`), or the donor's essential reactions for that carbon source (`essential`, precomputed by [`donor_essentiality_sweep.m`](scripts/gapfilling/curation/donor_essentiality_sweep.m)).

   &nbsp;

   b. **Scoring methods**\
   Gap-filling selection and curation score each model against the carbon source screen. Both parameters are set in [`gapfilling_selection.m`](scripts/gapfilling/selection/gapfilling_selection.m) and [`gapfilling_curation.m`](scripts/gapfilling/curation/gapfilling_curation.m).

   - `growthDataMode` sets the type of growth data: binary growth calls (`binary`) or continuous growth measurements (`continuous`).
   - `scoringMethod` sets the score. Binary data supports the Matthews correlation coefficient (`mcc`, default) and the area under the ROC curve (`auc`). Continuous data adds the Pearson correlation coefficient (`pearson`), normalized root-mean-square error (`nrmse`, `wnrmse`, `log-nrmse`), and symmetric mean absolute percentage error (`smape`).

   &nbsp;

5. **Outputs**\
   Generated models are located in [`data/models/`](data/models/), in the following structure:

   ```
   data/models/
   ├── draft/
   │   ├── carveme/             draft reconstructions from CarveMe (SBML)
   │   └── cobrapy/
   │       ├── mat/             draft reconstructions (.mat)
   │       └── xml/             draft reconstructions (SBML)
   ├── gapfilled/<method>/      gap-filled models, one folder per scoring method
   ├── curated/<method>/        curated models, one folder per scoring method
   └── final/
       ├── mat/                 final models (.mat)
       ├── xml/                 final models (SBML)
       └── reports/             MEMOTE quality reports (.html)
   ```

   Intermediate gap-filling and curation results are located in [`data/gapfilling/`](data/gapfilling/).

---

## Change history

This is version 2.0.1 of the pipeline. See [ChangeLog.md](ChangeLog.md) for a record of the methodological changes made since the previous version, [v1.0.1](https://github.com/VorholtLab/i-At-LSPHERE/releases/tag/v1.0.1) of [i-At-LSPHERE](https://github.com/VorholtLab/i-At-LSPHERE).

---

## References

1. Machado, D., Andrejev, S., Tramontano, M. & Patil, K. R. Fast automated reconstruction of genome-scale metabolic models for microbial species and communities. *Nucleic Acids Res.* **46**, 7542–7553 (2018). [https://doi.org/10.1093/nar/gky537](https://doi.org/10.1093/nar/gky537)
2. Vayena, E. *et al.* A workflow for annotating the knowledge gaps in metabolic reconstructions using known and hypothetical reactions. *Proc. Natl Acad. Sci. USA* **119**, e2211197119 (2022). [https://doi.org/10.1073/pnas.2211197119](https://doi.org/10.1073/pnas.2211197119)
3. Schäfer, M. *et al.* Metabolic interaction models recapitulate leaf microbiota ecology. *Science* **381**, eadf5121 (2023). [https://doi.org/10.1126/science.adf5121](https://doi.org/10.1126/science.adf5121)
4. Lieven, C. *et al.* MEMOTE for standardized genome-scale metabolic model testing. *Nat. Biotechnol.* **38**, 272–276 (2020). [https://doi.org/10.1038/s41587-020-0446-y](https://doi.org/10.1038/s41587-020-0446-y)
5. Bai, Y. *et al.* Functional overlap of the *Arabidopsis* leaf and root microbiota. *Nature* **528**, 364–369 (2015). [https://doi.org/10.1038/nature16192](https://doi.org/10.1038/nature16192)
