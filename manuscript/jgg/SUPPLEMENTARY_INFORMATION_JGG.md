# Supplementary Information

## Supplementary Results

### External-cohort detail

The animal-level external validation plots and continuous-effect comparisons are provided in Fig. S1. The 25 externally meta-replicated genes include p53-associated injury genes (*Cdkn1a*, *Gdf15*, *Phlda3*, *Zmat3* and *Bax*), remodeling genes (*Ecm1*, *Hapln1*, *Tinagl1*, *Serpine2* and *Lgals1*) and vascular genes (*Apln* and *Ddah1*). The pathway consensus links oxidative phosphorylation, cellular respiration, ATP metabolism and mitochondrial translation with reactive-oxygen-species–NFE2L2 signaling, mTORC1 activity, EMT-like signaling and wound healing. Figure S1 documents the study-specific effect dispersion around these shared directions.

### Exact animal-subset sensitivity in GSE243129

GSE243129 contained two animals in each age-by-oxygen stratum. All 81 combinations of non-empty normoxia and hyperoxia subsets, 16 one-pair-per-stratum configurations and eight leave-one-animal configurations were enumerated with fixed model predictions. scGen-style prediction ranked first in 67/81 configurations, 9/16 one-pair configurations and all eight leave-one-animal analyses. These analyses quantify influence within the deposited eight animals. Biological replication remains two animals per age-by-oxygen stratum (Fig. S2).

### Network-constrained perturbation priorities

Animal-balanced hyperoxia effects projected onto CollecTRI regulons yielded recurrent positive activities for Egr1, Srf, Smad4, Jund, Hif1a, Ets2, Sp1, Trp73, Trp63, Atf2, Trp53 and Jun in capillary subtypes (Fig. S3A). Across 329 shared transcription factors, the mean P14 capillary activity profile correlated with the independent hyperoxia contrast (Pearson r = 0.618; Spearman ρ = 0.640). Removing fitted single-factor contributions produced small changes in the 33-gene effect magnitude; the largest median reduction was 1.33% for Egr1. Several regulators jointly carry the response.

Trp53 activity rose in five of six age-by-capillary contrasts and was highest at P14. The independent hyperoxia contrast also increased Trp53 activity (ULM score 5.29, FDR = 4.70×10−6). Global p53 loss shifted the score negative and reversed the independent hyperoxia direction for 16 of 28 evaluable genes from the 33-gene set. Endothelial-specific p53 deletion retained a positive score and reversed six genes (Fig. S3B and S3C). These differences place p53 within a context-dependent stress network. Animal-balanced ligand–receptor analysis identified recurrent Gdf15–Tgfbr2, Apln–Aplnr, Tgfb2–Acvr1/Tgfbr1, Adm–Calcrl and Ccl4–Ackr2 routes. These expression-supported hypotheses connect the replicated injury program to capillary maintenance, TGF-β-responsive remodeling and vascular tone (Fig. S3D).

### Exploratory four-donor human projection

The original four-donor GSE275938 projection produced positive mouse-oriented scores in BPD-associated pulmonary hypertension compared with BPD alone. Equal-cell, equal-UMI and matched-random-set analyses retained that pattern. Diagnosis was aligned with donor recovery and sequencing depth, so the result remains a technically confounded observation. The donor-rich LungMAP atlas supplied the human test.

## Supplementary Figure legends

### Fig. S1. External mouse validation at animal, gene-set and pathway levels.

A: Locked 33-gene scores for room-air and hyperoxia animals in GSE243129, GSE209664 and the unbiased male-only GSE230672 contrast. B: External versus discovery log2 fold-changes for locked genes. C: Direction-concordant fractions for the 33- and 198-gene sets. D: Selected locked-pathway effects across the three designated external contrasts. RA, room air; HOX, hyperoxia; NES, normalized enrichment score.

### Fig. S2. Exact animal-subset sensitivity of the eight-animal GSE243129 benchmark.

Model ranks were recalculated over all non-empty animal subsets, one-pair-per-stratum configurations and leave-one-animal configurations with fixed predictions. Blue, green and orange bars denote 81 non-empty subsets, 16 one-pair-per-stratum configurations and eight leave-one-animal configurations, respectively. The analysis evaluates dependence on the eight deposited animals; biological replication remains two animals per age-by-oxygen stratum (eight animals in total).

### Fig. S3. Network-constrained perturbation priorities.

Animal-balanced effects were projected onto signed CollecTRI regulons and expression-supported OmniPath ligand–receptor routes. Single-factor counterfactual projections rank perturbation candidates. Panel D lists recurrent routes targeting capillary endothelium, including Gdf15–Tgfbr2, Apln–Aplnr and Tgfb2 receptor signaling. Gold bars denote routes also reported in the independent p53 CellChat analysis; purple bars denote the remaining recurrent routes. TF, transcription factor; ULM, univariate linear model.

### Fig. S4. Exploratory four-donor GSE275938 projection.

Mouse-oriented scores are shown across deposited endothelial subtypes and sampling sensitivities. Diagnosis was aligned with donor recovery and sequencing depth; this dataset is retained as a technically confounded observation and not as human population validation. EC, endothelial cells; aCap, aerocyte capillary; gCap, general capillary; BPD+PH, bronchopulmonary dysplasia with pulmonary hypertension; UMI, unique molecular identifier.
