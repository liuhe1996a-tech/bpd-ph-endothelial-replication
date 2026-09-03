# Cross-cohort single-cell genomics of neonatal hyperoxia reveals reproducible endothelial responses and cohort-dependent model rankings

## Abstract

Neonatal hyperoxia disrupts pulmonary microvascular development, yet the endothelial responses that persist across studies remain uncertain. We integrated animal-level discovery, internal replication and three external mouse cohorts, then evaluated perturbation prediction in two independently trained single-cell datasets and mouse-to-human concordance in a 21-donor LungMAP atlas. Twenty-five of 29 evaluable multi-subtype genes retained direction across all external mouse cohorts and passed signed meta-analysis FDR control. The replicated genes included *Cdkn1a*, *Gdf15*, *Phlda3*, *Zmat3*, *Bax*, *Ecm1* and *Apln*. Concordant pathways linked oxidative phosphorylation and ATP metabolism with reactive-oxygen-species–NFE2L2 signaling, mTORC1 activity and matrix remodeling. Among seven prediction methods, Sinkhorn optimal transport led in GSE151974 (mean capillary Spearman correlation, 0.534), while scGen-style prediction led in the eight-animal GSE243129 cohort (0.412). The LungMAP analysis showed near-zero locked-gene scores and no pathway-level directional replication; equivalence was supported within ±0.50 SD and remained unresolved within ±0.30 SD. Neonatal hyperoxia thus produced a reproducible mouse endothelial program, while predictive rankings changed by cohort and broad human concordance was weak.

**Keywords:** single-cell transcriptomics; perturbation prediction; functional genomics; biological replication; cross-species transfer; neonatal hyperoxia; endothelium

## Introduction

Bronchopulmonary dysplasia (BPD) develops during interrupted alveolar and pulmonary vascular maturation. Oxygen exposure, ventilation, inflammation and other perinatal insults alter the endothelial programs that support capillary growth and repair (Jobe and Bancalari, 2001; Jensen et al., 2019; Thébaud et al., 2019). Mouse studies have identified capillary subtype sensitivity, p53-associated injury responses and endothelial remodeling, while human BPD with pulmonary hypertension has been linked to loss of semaphorin signaling and FOXF1 activity (Hansmann et al., 2021; Hurskainen et al., 2021; Zanini et al., 2023; Shirazi et al., 2025). The signals that reproduce across independent hyperoxia cohorts, and their relationship to human disease, remain unresolved.

Single-cell data can localize these responses to endothelial subtypes and developmental stages. Their statistical resolution is still set by the number of animals or donors. Cells from one biological unit share exposure, tissue processing and biological history; reuse of those cells across feature selection, fitting and validation inflates apparent generalization (Squair et al., 2021). We addressed this problem at the study-design stage by assigning each dataset a fixed analytical role. Inference was performed at the animal or donor level.

Perturbation-prediction models offer a second view of the same biology. Linear shifts, latent-state arithmetic and optimal-transport maps estimate a hyperoxia response from normoxic cells, yet recent benchmarks report substantial dataset and endpoint dependence (Ahlmann-Eltze et al., 2025; Mesue Njume et al., 2026; Wei et al., 2026). We compared seven approaches with fold-specific features, animal-grouped validation, repeated stochastic fitting and shared biological resampling in two neonatal-hyperoxia cohorts.

The analysis combined purified endothelial RNA sequencing, animal-by-subtype single-cell pseudobulks, three unused mouse validation cohorts, p53 perturbation data and a donor-level LungMAP atlas (Fig. 1; Table 1; Table S1). We asked which endothelial genes and pathways survive cross-cohort testing, whether model rankings persist in an independently trained cohort, and whether the locked mouse program is evident in human BPD-associated endothelium. The principal biological signals are summarized in Table 1.

[Figure 1 near here]

[Table 1 near here]

## Results

### Neonatal hyperoxia induces a multi-subtype endothelial stress and remodeling program

GSE216046 supplied four air and four hyperoxia P14 mice. GSE151974 supplied 61,839 lung cells from 36 mice at P3, P7 and P14, including 10,969 endothelial cells. After filtering, 14,561 genes were modeled in purified P14 endothelium; 5,543 had FDR < 0.05, 1,029 also had observed |log2 fold-change|≥1, and 341 passed the formal effect-threshold test (Fig. 2A). In the P14 internal pseudobulk contrast, 174 of 179 significant capillary overlaps and 127 of 128 aerocyte-capillary overlaps agreed in direction. Across five endothelial subtypes, 198 genes met the internal replication rule and 33 reproduced in at least three subtypes (Fig. 2B and 2C). Of those 33, 21 passed the formal DESeq2 threshold, 15 passed an edgeR effect-threshold test in at least three subtypes, and 11 passed both strict criteria.

The 33-gene set resolved linked components of the endothelial response. *Cdkn1a*, *Gdf15*, *Phlda3*, *Zmat3*, *Bax* and *Rps27l* marked a p53-associated injury program. *Ecm1*, *Hapln1*, *Tinagl1*, *Abi3bp*, *Serpine2* and *Lgals1* connected hyperoxia to matrix and wound remodeling, while *Apln* and *Ddah1* retained capillary vascular relevance. The pathway analysis identified 162 directionally replicated terms: 95 Gene Ontology Biological Process, 19 Hallmark and 48 Reactome pathways. Respiratory electron transport, oxidative phosphorylation and complex-I assembly reproduced across all five endothelial subtypes. Reactive-oxygen-species, p53, apoptosis, mTORC1 and extracellular-matrix programs occupied narrower subtype distributions (Fig. 2D).

The broader 198-gene response was enriched in capillary cells in a post-selection, within-animal analysis. The age-adjusted capillary-minus-other oxygen effect was 0.129 signed-score units (95% CI 0.036–0.223; HC3 P = 0.0068; stratified-randomization P = 0.0063). The 33-gene contrast was smaller (0.027, 95% CI −0.089 to 0.144; HC3 P = 0.648; randomization P = 0.644). Capillary enrichment characterized the broad replicated response in this internal dataset; the compact multi-subtype set showed a more evenly distributed endothelial pattern.

[Figure 2 near here]

### Independent mouse cohorts retain p53-associated injury, bioenergetic and remodeling signals

Three previously unused mouse cohorts tested the locked evidence (Fig. S1). For the 198-gene union, 84 genes were evaluable in every external contrast; 68 followed the discovery direction in all three and 66 passed the signed meta-analysis FDR criterion. For the 33-gene set, 29 were evaluable in every contrast, 25 agreed in all three and all 25 passed the signed meta-analysis criterion (Fig. 3A–3C). The retained genes included the injury-response genes *Cdkn1a*, *Gdf15*, *Phlda3*, *Zmat3* and *Bax*; the remodeling genes *Ecm1*, *Hapln1*, *Tinagl1*, *Serpine2* and *Lgals1*; and the vascular genes *Apln* and *Ddah1*. *Cyp26b1*, *Sox11*, *Abi3bp*, *Hoxa5*, *Tox3* and *Zfp36l2* reproduced as downregulated signals.

The pathway results connected this gene-level response to mitochondrial bioenergetics, oxidative stress and tissue remodeling. With one designated comparison per study, 71 of 88 evaluable Gene Ontology, 13 of 19 Hallmark and 38 of 47 Reactome terms were directionally concordant in all three external cohorts; 22, 4 and 8, respectively, also had FDR < 0.05 in every cohort (Fig. 3D). Oxidative phosphorylation, cellular respiration, ATP metabolism and mitochondrial translation reproduced alongside reactive-oxygen-species and NFE2L2 signaling, mTORC1 activity, epithelial–mesenchymal-transition-like signaling and wound healing. The independent cohorts retained a coordinated mouse endothelial program spanning energy metabolism, cellular stress and structural remodeling (Fig. 3A–3D; Table S2).

[Figure 3 near here]

### Prediction point-estimate rankings differ across independently trained cohorts

Seven methods were evaluated on animal-level hyperoxia effects at a completely held-out age. In the primary GSE151974 cohort, mean Spearman correlation across six capillary tasks was 0.251 for gene-space shift, 0.514 for PCA shift, 0.491 for VAE, 0.440 for conditional VAE, 0.464 for scGen-style, 0.470 for CPA-style and 0.534 for Sinkhorn optimal transport (Fig. 4A). Across 500 joint-animal bootstrap replicates, Sinkhorn exceeded PCA by 0.0162 (95% percentile interval 0.0099 to 0.0231). VAE, conditional VAE, scGen-style and CPA-style underperformed PCA by −0.0197 (−0.0306 to −0.0066), −0.0618 (−0.0820 to −0.0425), −0.0413 (−0.0564 to −0.0237) and −0.0370 (−0.0521 to −0.0195), respectively (Fig. 4B). Seed-level standard deviations ranged from 0.0040 for Sinkhorn to 0.0162 for CPA-style.

Point-estimate ordering differed in the independently trained GSE243129 cohort, which contained only two animals per age-by-oxygen group. Across four held-out age-by-capillary tasks, mean Spearman correlations were 0.043 for gene-space shift, 0.364 for PCA, 0.405 for VAE, 0.409 for conditional VAE, 0.412 for scGen-style, 0.403 for CPA-style and 0.384 for Sinkhorn (Fig. 4C). Under this resampling scheme, all five stochastic methods had higher point estimates than PCA: median paired differences were 0.0415 for VAE (0.0225 to 0.0666), 0.0451 for conditional VAE (0.0136 to 0.0752), 0.0473 for scGen-style (0.0251 to 0.0753), 0.0382 for CPA-style (0.0073 to 0.0570) and 0.0170 for Sinkhorn (0.0100 to 0.0246). After draw order was removed, the 500 bootstrap replicates represented 77 of 81 possible unordered resampling configurations. In the exact non-empty-subset analysis, scGen-style ranked first in 67/81 configurations and in all eight leave-one-animal analyses; it ranked first in 9/16 one-pair-per-stratum configurations (Fig. S2). These sensitivity results describe ordering stability within the eight deposited animals. Biological replication remained two animals per age-by-oxygen group. The primary-cohort point-estimate leader differed from that in this small secondary cohort (Fig. 4A–4D; Fig. S2; Table S3).

### Evaluation endpoints alter the apparent advantage of model classes

Endpoint choice further changed the apparent method advantage. On the 198-gene–HVG intersection in GSE151974, VAE and PCA were nearly identical (median difference −0.0003, −0.0097 to 0.0103), and Sinkhorn exceeded PCA by 0.0032 (−0.0009 to 0.0088). The 33-gene–HVG intersection contained only 12–13 genes per fold and produced unstable exploratory estimates. Independent GSE266988 P14 calibration also showed method- and subtype-dependent performance: on 772 shared primary genes, mean Spearman correlations across general and aerocyte capillary cells were 0.397 for gene-space, 0.436 for PCA, 0.417 for VAE, 0.331 for conditional VAE, 0.375 for scGen-style, 0.396 for CPA-style and 0.451 for Sinkhorn. Model class alone did not determine performance across cohorts and endpoints.

[Figure 4 near here]

### Human BPD-associated endothelium shows little concordance with the locked mouse program



The independent LungMAP atlas provided a donor-level test with 13 primary BPD and 8 primary normal donors. All BPD donors were born preterm and all controls were born at term, so the fitted coefficient represents the combined prematurity/BPD state. In all endothelium, the adjusted strict 25-gene score difference was 0.037 (95% CI −0.321 to 0.395; HC3 P = 0.830; permutation P = 0.845), and the extended 33-gene difference was 0.010 (95% CI −0.323 to 0.343; P = 0.951; permutation P = 0.951; Fig. 5A). The strict-score estimate was also near zero in combined capillary cells (0.014, 95% CI −0.339 to 0.368), CAP1 (0.073, −0.288 to 0.435) and CAP2 (0.085, −0.214 to 0.385). None of the eight primary signature tests approached the multiplicity-adjusted criterion. For the strict score, equivalence was supported at ±0.50 standard deviations in all endothelium (TOST P = 0.0072) and combined capillary cells (P = 0.0050); the ±0.30-standard-deviation margins remained unresolved (P = 0.070 and 0.053, respectively). The approximate 80%-power minimum detectable effects were 0.476 and 0.470 standard deviations. These estimates exclude a broad average effect of at least 0.5 standard deviations under the combined state contrast and leave smaller effects unresolved.

Gene-level tests did not support directional transfer. In all endothelium, 13 of 25 strict genes followed the mouse direction (one-sided binomial P = 0.500), and the mouse discovery effects were uncorrelated with adjusted human effects (Spearman ρ = −0.011, P = 0.959; Fig. 5B). Expression- and detection-matched random-set tests were non-significant in every population and signature (empirical FDR = 0.727). Leave-one-donor-out estimates remained small; stage-specific models had minimum unadjusted P = 0.073 and minimum FDR = 0.989. Endothelium-versus-epithelial, mesenchymal or immune interactions also failed FDR control.

Sampling sensitivities remained centered near zero (Fig. 5C). Equal-UMI intervals crossed zero for every population and signature. Equal-cell sampling was centered near zero or negative; with only 15 cells per donor, the CAP1 median was −0.115 for the strict score and −0.105 for the extended score. Of 162 locked mouse pathways, 156 had same-name human counterparts. None met the joint requirement of the mouse direction and human FDR < 0.05 in all endothelium, combined capillary cells, CAP1 or CAP2 (Fig. 5D). Only 42 of 156 pathways in all endothelium followed the mouse direction, with similarly low fractions in the capillary analyses. The larger human atlas showed little broad directional concordance with the mouse hyperoxia program; effects near 0.3 standard deviations remain within the statistical uncertainty (Fig. 5A–5D; Table S4).

[Figure 5 near here]

## Discussion

The strongest biological result is a mouse endothelial program that survived internal replication and three independent cohort tests. Its gene-level structure connected p53-associated injury (*Cdkn1a*, *Gdf15*, *Phlda3*, *Zmat3* and *Bax*) with extracellular-matrix and wound-remodeling genes (*Ecm1*, *Hapln1*, *Tinagl1*, *Serpine2* and *Lgals1*) and capillary vascular signals such as *Apln* and *Ddah1*. The pathway analysis placed these genes within a broader response involving oxidative phosphorylation, ATP metabolism, mitochondrial translation, reactive-oxygen-species–NFE2L2 signaling, mTORC1 activity and matrix remodeling. This combination indicates that neonatal oxygen injury recruits metabolic adaptation and structural repair alongside canonical stress signaling. The capillary enrichment of the 198-gene response further links the program to the endothelial compartment most directly involved in alveolar gas-exchange development.

The recurrent regulatory and communication signals define specific experiments. *Trp53* activity increased in five of six age-by-capillary tasks and in an independent external contrast. Global p53 loss reversed more locked responses than endothelial-specific deletion, placing p53 within a distributed, context-dependent injury network. The recurrent Gdf15–Tgfbr2 route connects the cross-cohort injury marker *Gdf15* to TGF-β-responsive remodeling. Apln–Aplnr retains a capillary homeostatic axis, while Tgfb2–Acvr1/Tgfbr1 and Adm–Calcrl extend the candidate network toward vascular remodeling and tone. These expression-supported routes can now be tested with endothelial perturbation, receptor blockade and vascular phenotyping.

The prediction benchmark adds a distinct functional-genomics result. Sinkhorn optimal transport produced the highest point estimate in GSE151974, while scGen-style prediction led in GSE243129. Endpoint restriction to the internally selected disease genes further compressed or changed method differences. All models used fold-specific features, animal-grouped validation, ten stochastic seeds where applicable and shared animal resampling. The ordering changed with the cohort and evaluation endpoint under the same protocol. The eight-animal secondary cohort supports a hypothesis-generating comparison; its two animals per age-by-oxygen group set the biological resolution of that result.

The LungMAP analysis placed a clear boundary on cross-species interpretation. The locked mouse genes showed a near-zero average score, no gene-level correlation and no pathway meeting the joint direction-and-FDR rule in the primary human populations. Equivalence testing excluded an average effect of at least 0.5 standard deviations and left effects near 0.3 standard deviations unresolved. Human BPD develops across prematurity, ventilation, inflammation, infection, nutrition and treatment, whereas the mouse cohorts model a defined oxygen exposure. In the LungMAP atlas, gestational age and BPD status were completely aligned, making the human coefficient a combined prematurity/BPD-state estimate. The four-donor GSE275938 analysis reinforces this caution: a visually positive score coincided with differential cell recovery and sequencing depth (Fig. S4).

Several limitations define the scope of the conclusions. All data were public and no new wet-laboratory experiment was performed. Bulk and single-cell platforms differ, deposited cell labels were accepted after audit, and small external cohorts limit interaction tests. The prediction benchmark covers two related neonatal-hyperoxia cohorts. The protocol-adapted scGen and CPA implementations reproduce their defining latent operations under a common split and tuning budget; they do not reproduce every feature of the official software. Sinkhorn is an entropy-regularized transport baseline and differs from the parameterized neural CellOT model. Foundation models were outside the benchmark. Hyperparameters were fixed to protect held-out outcomes from tuning. The LungMAP comparison could not separate prematurity from BPD, and exact-name pathway mapping left six locked pathways untestable.

Neonatal hyperoxia produces a reproducible mouse endothelial program that couples mitochondrial bioenergetics, oxidative stress, p53-associated injury and matrix remodeling. The program prioritizes GDF15–TGFBR2, APLN–APLNR and TGF-β receptor signaling for functional testing. Prediction performance depends on the cohort and evaluation target, and the 21-donor human atlas shows little evidence for broad directional transfer. These findings convert a heterogeneous set of public datasets into a focused molecular model with explicit experimental priorities and defined limits of human translation.

## Materials and methods

### Study design and evidence locks

This secondary study analyzed public, de-identified transcriptomic data. GSE216046 was assigned to mouse discovery (Long et al., 2022; NCBI, 2026c). GSE151974 was assigned to internal animal-by-subtype replication, feature locking and the primary perturbation-prediction benchmark (Hurskainen et al., 2021; NCBI, 2026a). GSE243129, GSE209664 and GSE230672 were not used to select biological signatures and served as external mouse validation cohorts (Gong et al., 2020; NCBI, 2026b; NCBI, 2026d; NCBI, 2026e); GSE243129 additionally supplied a second, independently trained perturbation-prediction benchmark. Published GSE266988 capillary and p53-perturbation contrasts were used for external calibration and mechanistic triangulation (Vila Ellis et al., 2025; NCBI, 2026f). GSE275938 was retained as a four-donor exploratory BPD-associated pulmonary hypertension projection (Shirazi et al., 2025; NCBI, 2026g). A versioned LungMAP/CELLxGENE BPD atlas was the independent human transfer test (LungMAP Consortium, 2026). The statistical unit was the sample, animal or donor, never the cell.

The 198-gene replicated union, 33-gene multi-subtype set and 162 mouse-replicated pathways were frozen before external mouse outcomes and all human outcomes were summarized. Membership remained fixed through external testing. Within the 33-gene set, 25 genes that subsequently met the external mouse meta-replication rule formed a stricter cross-mouse signature. This 25-gene subset was fixed before analysis of the LungMAP atlas.

### Mouse discovery and internal single-cell replication

For GSE216046, genes were retained at counts per million of at least 1 in at least three libraries and analyzed with DESeq2 (Love et al., 2014). The design contrasted hyperoxia with air. Benjamini–Hochberg false-discovery rate (FDR) control was applied (Benjamini and Hochberg, 1995). FDR < 0.05 with an observed absolute log2 fold-change of at least 1 defined the descriptive discovery screen; a formal DESeq2 threshold test against |log2 fold-change|>1 was a stricter sensitivity analysis. Variance-stabilized counts were used only for principal-component analysis.

Raw GSE151974 UMI counts were summed within animal-by-endothelial-subtype groups. edgeR filtering, TMM normalization, robust dispersion estimation and quasi-likelihood fitting were used (Lun et al., 2016; Chen et al., 2025). The primary model included categorical age, oxygen and their interaction; the postnatal day 14 (P14) contrast combined the oxygen main effect and P14 interaction. P14-only, all-age, minimum-cell and subtype-exclusion specifications were sensitivity analyses. A gene was internally replicated within a subtype if it was measured in both datasets, passed the discovery screen, had GSE151974 FDR < 0.05 and had the same direction. The union across five endothelial subtypes defined the 198-gene set; genes reproduced in at least three subtypes defined the 33-gene set.

Genes were ranked by signed model statistics for gene-set enrichment analysis of mouse MSigDB 2026.1 Hallmark, Reactome and Gene Ontology Biological Process collections (Subramanian et al., 2005; Liberzon et al., 2015; Korotkevich et al., 2019; Castanza et al., 2023). Pathway replication required FDR < 0.05 and the same normalized-enrichment-score direction in both discovery and internal replication. Composition was analyzed per animal with exact two-sided oxygen-label permutations and FDR correction. Recovered-cell fractions were treated as composition measurements. Capillary enrichment was tested with signed gene scores calculated in animal-by-subtype pseudobulks. The within-animal contrast was the mean of general and aerocyte capillary scores minus the mean of arterial, venous and lymphatic scores. An age-adjusted oxygen coefficient with HC3 standard errors and 100,000 oxygen-label randomizations within age strata tested whether hyperoxia produced a larger signed response in capillary than in other endothelial subtypes. The 198-gene and 33-gene endpoints were tested separately.

### Independent mouse validation and signed evidence synthesis

External datasets were processed without refitting the locked features. GSE243129 analyses retained wild-type animals and used animal-by-endothelial-compartment pseudobulks at P7 and P14. Its contribution to the signed gene meta-analysis was the P14 general-capillary (gCap) hyperoxia-versus-normoxia contrast from four animals; P7 and aerocyte-capillary (Cap-a) effects were not pooled or averaged. GSE209664 models retained all 12 animals and adjusted oxygen contrasts for sex; oxygen labels were permuted within sex for exact signature and composition tests. GSE230672 contained two hyperoxia females, two hyperoxia males and two room-air males. Because exposure and sex were inseparable among females, the primary contrast was the unbiased two-versus-two male comparison. The full imbalanced dataset was descriptive only.

For each gene, a two-sided external P value was converted to a signed z score using the direction of its external log2 fold-change. Signed z scores were combined with square-root sample-size weights. Combined P values were two-sided and adjusted across locked genes. External meta-replication required the external effect to follow the discovery direction in every evaluable cohort and the signed meta-analysis to retain that direction at FDR < 0.05. Pathways were evaluated separately within the three locked MSigDB collections. Each external study contributed exactly one designated effect per gene and one pathway contrast: P14 gCap hyperoxia versus normoxia for GSE243129, P14 capillary hyperoxia versus room air adjusted for sex for GSE209664, and P14 male sorted endothelium hyperoxia versus room air for GSE230672. This designation prevented ages or endothelial compartments from the same study from contributing multiple effects. Pathway consensus required directional concordance in all three studies.

### Perturbation-prediction inputs, models and held-out evaluation

The deposited GSE151974 matrix was restricted to 10,969 endothelial cells. Counts were library-size normalized to 10,000 per cell and log1p transformed. Three outer folds held out P3, P7 or P14. Every cell from the held-out age, including normoxia and hyperoxia cells, was excluded from feature selection, scaling, principal-component fitting, neural-network fitting and epoch selection. Within each outer-training fold, the model matrix contained exactly 1,800 highly variable genes selected among genes detected in at least max(20, 0.5% of training cells). No member of the internally selected 198- or 33-gene sets was added unless it independently satisfied the fold-specific highly variable-gene rule. The full fold-specific highly variable-gene set was the primary endpoint. Intersections with the 198- and 33-gene sets were reported only as post-selection, disease-informed sensitivity analyses because their membership had been selected using GSE151974. Held-out normoxia cells were used only after fitting as the starting state for prediction.

Seven predictive approaches were compared. The gene-space baseline added the animal-balanced hyperoxia-minus-normoxia response learned from the other ages to held-out normoxia animal profiles. The principal-component baseline learned the same shift in 32 components. The VAE used 256- and 128-unit GELU encoder layers, a 32-dimensional latent space and a symmetric decoder. A conditional VAE added scaled age, subtype and oxygen indicators to its encoder and decoder. A protocol-adapted scGen model combined the same VAE backbone with animal-, age- and subtype-balanced latent vector arithmetic (Lotfollahi et al., 2019). A protocol-adapted compositional perturbation autoencoder (CPA) learned a basal representation and additive oxygen and context embeddings, with gradient-reversal adversaries discouraging residual perturbation and context information (Lotfollahi et al., 2023). Both implementations used the shared architecture, split and tuning budget and are described as scGen-style and CPA-style. The optimal-transport baseline fitted entropy-regularized Sinkhorn couplings between animal-balanced normoxia and hyperoxia cells in outer-fold-trained principal-component space, then used a 15-nearest-neighbor out-of-sample displacement map for held-out normoxia cells. This explicit transport baseline evaluates an optimal-transport prior and differs from CellOT's parameterized neural map (Bunne et al., 2023).

All seven methods received the same outer-training cells and the same 1,800 fold-specific genes. Neural models used AdamW with learning rate 0.001, weight decay 10−5, dropout 0.1, KL weight 10−4 and a maximum of 60 epochs. For epoch selection, one complete animal was assigned to validation from each outer-training age-by-oxygen stratum; no animal contributed cells to both fitting and validation. The selected epoch was then used to refit the model on all outer-training animals. VAE, conditional VAE, scGen-style, CPA-style and Sinkhorn predictions were repeated for ten fixed analysis seeds (20260817–20260826); seeds controlled initialization or balanced cell subsampling as applicable. A no-change prediction was retained as an audit baseline.

Predicted cells were aggregated to animal means before evaluation. Pearson and Spearman correlations, root-mean-square error and direction accuracy were calculated for the fold-specific primary genes and the two disease-set intersections. For each held-out age, five hundred bootstrap replicates drew normoxia and hyperoxia animal identifiers once and applied those same identifiers to both capillary subtypes and every model; each predicted profile remained paired to its starting normoxia animal. Stochastic predictions were averaged across the ten seeds within each biological bootstrap replicate; seed-to-seed variation was summarized separately. Paired model-minus-baseline differences were then averaged across six capillary tasks (two capillary subtypes across three ages), and percentile intervals over the 500 joint biological resamples were reported.

GSE243129 was analyzed as a second prediction cohort with cohort-specific model fitting. Its eight wild-type animals comprised two animals in each age-by-oxygen group at P7 and P14 and supplied four held-out age-by-capillary tasks. Each fold selected exactly 1,800 HVGs from the other age, grouped validation by animal, used the same architectures, fixed tuning budget, ten seeds and 500 joint-animal bootstrap replicates, and imported no GSE151974 features, scaling parameters or model weights. The bootstrap quantified resampling variability from the eight deposited animals. We additionally enumerated all 81 combinations of non-empty normoxia and hyperoxia animal subsets across the four age-by-oxygen strata, 16 configurations retaining one animal per stratum, and eight leave-one-animal configurations; each age-specific selection was applied jointly to Cap and Cap-a, model predictions remained fixed, and ranks were calculated after averaging stochastic seeds and the four held-out tasks. This analysis replicated the evaluation protocol with models trained independently in GSE243129. Independent GSE266988 P14 capillary effects were intersected by gene and scored without model selection or refitting as an additional calibration analysis.

### Network-constrained prioritization

For each age and endothelial subtype, log-normalized expression was first averaged within animal and then across animals to form hyperoxia-minus-normoxia effect vectors. Mouse CollecTRI signed transcription-factor–target interactions were analyzed with the univariate linear model in decoupler, requiring at least ten targets (Badia-i-Mompel et al., 2022; Müller-Dott et al., 2023). Ordinary least-squares coefficients were fitted against signed regulon weights. Each fitted transcription-factor contribution was removed from the effect vector, and the proportional reduction in root-mean-square magnitude of the locked gene sets was recorded. These counterfactual projections rank perturbation candidates for experimental testing.

All 61,839 GSE151974 lung cells were also mapped to capillary, other endothelial, epithelial, fibroblast, mural, myeloid, lymphoid, mast or mesothelial compartments. OmniPath ligand–receptor genes were aggregated to animal-by-compartment pseudobulks (Türei et al., 2016). Expression, cell-count and route-coverage thresholds were applied before exact six-versus-six oxygen-label tests. Recurrence required a positive effect at two or more ages and nominal exact P ≤ 0.05 at one or more ages. Published GSE266988 CellChat pairs were used only as external annotations.

### Exploratory four-donor human projection

Mouse symbols were reconciled with NCBI Gene, and one-to-one mouse–human orthologues were fixed before reading human outcomes (Maglott et al., 2011). Raw GSE275938 counts were summed within donor-by-subtype groups. TMM log2 counts per million were oriented by the fixed mouse direction, standardized across donors and averaged with equal gene weights. The dataset contained two BPD and two BPD-associated pulmonary hypertension donors. Equal-cell, equal-total-UMI and expression-matched random-set analyses assessed technical sensitivity. Diagnosis was confounded with cell recovery and sequencing depth, so GSE275938 was retained as a technical sensitivity analysis. LungMAP supplied the donor-level human test.

### Independent LungMAP human transfer test

The versioned LungMAP H5AD object contained 271,381 nuclei from 24 donors and 35,477 feature rows (LungMAP Consortium, 2026). Raw integer counts were read from raw.X. Sixteen duplicate feature rows representing eight symbols were collapsed by summation, leaving 35,469 unique symbols; none involved a locked signature gene. The primary analysis followed the deposited is_primary_data flag and included 13 BPD and 8 normal donors, comprising 255,204 nuclei and 12,758 endothelial nuclei. Three additional normal donors were reserved for sensitivity analysis; the complete object contained 14,115 endothelial nuclei. We evaluated all endothelium, combined CAP1 and CAP2, CAP1, and CAP2. Epithelial, mesenchymal and immune lineages were negative-control comparators.

Counts were summed within donor and population and transformed to log2 counts per million with a 0.5 count offset. All 33 one-to-one orthologues were present. The strict 25-gene and extended 33-gene scores used mouse discovery directions, donor-standardized gene expression and equal gene weights. Primary ordinary least-squares models included disease state, postnatal age and sex; HC3 heteroscedasticity-consistent standard errors were used. Gestational age was omitted because all 13 primary BPD donors were born at 23–28 weeks and all eight primary controls were born at 40 weeks. The disease coefficient estimates the combined prematurity/BPD state. Freedman–Lane residual permutation provided 10,000 two-sided permutation tests. As a final post hoc sensitivity analysis, two one-sided tests evaluated equivalence at margins of ±0.30 and ±0.50 standard deviations; normal-approximation minimum detectable effects at 80% power were reported to distinguish absence of evidence from evidence against a large effect.

Sensitivity analyses included leave-one-donor-out refitting, all 24 donors, a minimum 20-cell pseudobulk rule, a postnatal-age-overlap subset, author-deposited BPD stages, 200 equal-cell samples and 200 equal-UMI rarefactions. Ten thousand random gene sets matched each signature gene on mean expression and detection within a filtered population-specific universe. Endothelium-versus-lineage contrasts tested whether any disease-oriented shift was endothelial-selective.

The 162 locked mouse pathways were mapped by exact name to human MSigDB 2026.1; 156 had same-name human counterparts. Donor-level gene statistics were the adjusted disease coefficient divided by its HC3 standard error. Preranked enrichment used 2,000 permutations. Human pathway replication required both the locked mouse direction and FDR < 0.05. Human outcomes did not alter the locked pathway list.

### Statistics and reproducibility

Tests were two-sided unless direction was part of a defined replication or empirical-null rule. FDR < 0.05 was the multiplicity-adjusted criterion. Exact resampling P values used the plus-one correction. No method was used to predetermine sample size, and investigators were not blinded because public group labels were required. Complete derived tables, source data, random seeds, input checksums, an executed audit notebook, a technical report, environment manifests and code accompany the study.

## Ethics statement

This secondary study used only public, de-identified transcriptomic data and involved no new recruitment, intervention, biospecimen collection or animal experiment. The original studies reported their applicable approvals. No additional ethics approval or consent was required.

## Data availability

Mouse and human GEO datasets are available under GSE216046, GSE151974, GSE243129, GSE209664, GSE230672, GSE266988 and GSE275938 (NCBI, 2026a; NCBI, 2026b; NCBI, 2026c; NCBI, 2026d; NCBI, 2026e; NCBI, 2026f; NCBI, 2026g). The LungMAP atlas is available from CELLxGENE collection 3a5dbf8a-9b3e-4309-b4c5-d8a024f83734, dataset version 83bbeaaf-f5e0-42ac-be8a-dbc4e0c0d433 (LungMAP Consortium, 2026). Complete derived tables and numerical source data are contained in the fixed release assets associated with GitHub release v1.0.0-jgg and the companion Zenodo record 10.5281/zenodo.22263929.

## Code availability

The complete versioned workflow, fixed seeds, environment manifests, input checksums, executed audit notebook and clean-reproduction report are available at https://github.com/liuhe1996a-tech/bpd-ph-endothelial-replication (release v1.0.0-jgg) and archived in Zenodo 10.5281/zenodo.22263929.

## CRediT authorship contribution statement

H.L. contributed conceptualization, methodology, software, formal analysis, data curation, visualization and writing of the original draft. X.F. contributed validation, investigation, data verification and review and editing. L.S. contributed methodology, supervision, clinical interpretation and review and editing. M.Y. contributed conceptualization, supervision, project administration, clinical interpretation and review and editing. All authors reviewed and approved the manuscript and accept accountability for the work.

## Conflict of interest

The authors declare no competing interests.

## Acknowledgments

The authors thank the investigators, participants and animal-research teams who generated and shared the public datasets.

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

## Declaration of generative AI and AI-assisted technologies in the manuscript preparation process

During preparation of this work, the authors used OpenAI Codex to assist with code drafting, consistency checking, document organization and language editing. The authors reviewed and edited all outputs, verified the source data and numerical results, and take full responsibility for the content of the article. No generative AI system is listed as an author, and no generative AI was used to create or alter scientific images.

## References

Ahlmann-Eltze, C., Huber, W., Anders, S., 2025. Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. Nat. Methods 22, 1657–1661. https://doi.org/10.1038/s41592-025-02772-6.

Badia-i-Mompel, P., Vélez Santiago, J., Braunger, J., Geiss, C., Dimitrov, D., Müller-Dott, S., Taus, P., Dugourd, A., Holland, C.H., Ramirez Flores, R.O., et al., 2022. DecoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinform. Adv. 2, vbac016. https://doi.org/10.1093/bioadv/vbac016.

Benjamini, Y., Hochberg, Y., 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J. R. Stat. Soc. Ser. B 57, 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x.

Bunne, C., Stark, S.G., Gut, G., del Castillo, J.S., Levesque, M., Lehmann, K.V., Pelkmans, L., Krause, A., Rätsch, G., 2023. Learning single-cell perturbation responses using neural optimal transport. Nat. Methods 20, 1759–1768. https://doi.org/10.1038/s41592-023-01969-x.

Castanza, A.S., Recla, J.M., Eby, D., Thorvaldsdóttir, H., Bult, C.J., Mesirov, J.P., 2023. Extending support for mouse data in the molecular signatures database (MSigDB). Nat. Methods 20, 1619–1620. https://doi.org/10.1038/s41592-023-02014-7.

Chen, Y., Chen, L., Lun, A.T.L., Baldoni, P.L., Smyth, G.K., 2025. EdgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets. Nucleic Acids Res. 53, gkaf018. https://doi.org/10.1093/nar/gkaf018.

Gong, J., Feng, Z., Peterson, A.L., Carr, J.F., Vang, A., Braza, J., Choudhary, G., Dennery, P.A., Yao, H., 2020. Endothelial to mesenchymal transition during neonatal hyperoxia‐induced pulmonary hypertension. J. Pathol. 252, 411–422. https://doi.org/10.1002/path.5534.

Hansmann, G., Sallmon, H., Roehr, C.C., Kourembanas, S., Austin, E.D., Koestenberger, M., 2021. Pulmonary hypertension in bronchopulmonary dysplasia. Pediatr. Res. 89, 446–455. https://doi.org/10.1038/s41390-020-0993-4.

Hurskainen, M., Mižíková, I., Cook, D.P., Andersson, N., Cyr-Depauw, C., Lesage, F., Helle, E., Renesme, L., Jankov, R.P., Heikinheimo, M., et al., 2021. Single cell transcriptomic analysis of murine lung development on hyperoxia-induced damage. Nat. Commun. 12, 1565. https://doi.org/10.1038/s41467-021-21865-2.

Jensen, E.A., Dysart, K., Gantz, M.G., McDonald, S., Bamat, N.A., Keszler, M., Kirpalani, H., Laughon, M.M., Poindexter, B.B., Duncan, A.F., et al., 2019. The diagnosis of bronchopulmonary dysplasia in very preterm infants. An evidence-based approach. Am. J. Respir. Crit. Care Med. 200, 751–759. https://doi.org/10.1164/rccm.201812-2348OC.

Jobe, A.H., Bancalari, E., 2001. Bronchopulmonary dysplasia. Am. J. Respir. Crit. Care Med. 163, 1723–1729. https://doi.org/10.1164/ajrccm.163.7.2011060.

Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M.N., Sergushichev, A., 2019. Fast gene set enrichment analysis. bioRxiv 060012. https://doi.org/10.1101/060012.

Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J.P., Tamayo, P., 2015. The molecular signatures database hallmark gene set collection. Cell Syst. 1, 417–425. https://doi.org/10.1016/j.cels.2015.12.004.

Long, Y., Chen, H., Deng, J., Ning, J., Yang, P., Qiao, L., Cao, Z., 2022. Deficiency of endothelial FGFR1 alleviates hyperoxia-induced bronchopulmonary dysplasia in neonatal mice. Front. Pharmacol. 13, 1039103. https://doi.org/10.3389/fphar.2022.1039103.

Lotfollahi, M., Wolf, F.A., Theis, F.J., 2019. ScGen predicts single-cell perturbation responses. Nat. Methods 16, 715–721. https://doi.org/10.1038/s41592-019-0494-8.

Lotfollahi, M., Klimovskaia Susmelj, A., De Donno, C., Hetzel, L., Ji, Y., Ibarra, I.L., Srivatsan, S.R., Naghipourfar, M., Daza, R.M., Martin, B., et al., 2023. Predicting cellular responses to complex perturbations in high‐throughput screens. Mol. Syst. Biol. 19, e11517. https://doi.org/10.15252/msb.202211517.

Love, M.I., Huber, W., Anders, S., 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 15, 550. https://doi.org/10.1186/s13059-014-0550-8.

Lun, A.T.L., Chen, Y., Smyth, G.K., 2016. It’s DE-licious: a recipe for differential expression analyses of RNA-seq experiments using Quasi-Likelihood methods in edgeR. Methods Mol. Biol. 1418, 391–416. https://doi.org/10.1007/978-1-4939-3578-9_19.

LungMAP Consortium, 2026. Bronchopulmonary dysplasia single-nucleus RNA-sequencing collection. CELLxGENE collection 3a5dbf8a-9b3e-4309-b4c5-d8a024f83734; dataset version 83bbeaaf-f5e0-42ac-be8a-dbc4e0c0d433. https://cellxgene.cziscience.com/collections/3a5dbf8a-9b3e-4309-b4c5-d8a024f83734 (accessed 19 August 2026).

Maglott, D., Ostell, J., Pruitt, K.D., Tatusova, T., 2011. Entrez gene: gene-centered information at NCBI. Nucleic Acids Res. 39, D52–D57. https://doi.org/10.1093/nar/gkq1237.

Mesue Njume, C., Petracci, I., Bellini, S., Goljanek-Whysall, K., Quinlan, L.R., Fiszer, A., Borroni, B., Ghidoni, R., Kumbasar, A., Cakmak, A., 2026. When complexity does not pay: benchmarking deep learning and ensemble methods for biomarker discovery. Brief. Bioinform. 27, bbag211. https://doi.org/10.1093/bib/bbag211.

Müller-Dott, S., Tsirvouli, E., Vazquez, M., Ramirez Flores, R.O., Badia-i-Mompel, P., Fallegger, R., Türei, D., Lægreid, A., Saez-Rodriguez, J., 2023. Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities. Nucleic Acids Res. 51, 10934–10949. https://doi.org/10.1093/nar/gkad841.

National Center for Biotechnology Information, 2026a. Gene Expression Omnibus GSE151974. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151974 (accessed 19 August 2026).

National Center for Biotechnology Information, 2026b. Gene Expression Omnibus GSE209664. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE209664 (accessed 19 August 2026).

National Center for Biotechnology Information, 2026c. Gene Expression Omnibus GSE216046. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216046 (accessed 19 August 2026).

National Center for Biotechnology Information, 2026d. Gene Expression Omnibus GSE230672. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230672 (accessed 19 August 2026).

National Center for Biotechnology Information, 2026e. Gene Expression Omnibus GSE243129. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243129 (accessed 19 August 2026).

National Center for Biotechnology Information, 2026f. Gene Expression Omnibus GSE266988. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266988 (accessed 19 August 2026).

National Center for Biotechnology Information, 2026g. Gene Expression Omnibus GSE275938. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275938 (accessed 19 August 2026).

Shirazi, S.P., Negretti, N.M., Jetter, C.S., Sharkey, A.L., Garg, S., Kapp, M.E., Wilkins, D., Fortier, G., Mallapragada, S., Banovich, N.E., et al., 2025. Bronchopulmonary dysplasia with pulmonary hypertension associates with semaphorin signaling loss and functionally decreased FOXF1 expression. Nat. Commun. 16, 5004. https://doi.org/10.1038/s41467-025-60371-7.

Squair, J.W., Gautier, M., Kathe, C., Anderson, M.A., James, N.D., Hutson, T.H., Hudelle, R., Qaiser, T., Matson, K.J.E., Barraud, Q., et al., 2021. Confronting false discoveries in single-cell differential expression. Nat. Commun. 12, 5692. https://doi.org/10.1038/s41467-021-25960-2.

Subramanian, A., Tamayo, P., Mootha, V.K., Mukherjee, S., Ebert, B.L., Gillette, M.A., Paulovich, A., Pomeroy, S.L., Golub, T.R., Lander, E.S., et al., 2005. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc. Natl. Acad. Sci. U. S. A. 102, 15545–15550. https://doi.org/10.1073/pnas.0506580102.

Thébaud, B., Goss, K.N., Laughon, M., Whitsett, J.A., Abman, S.H., Steinhorn, R.H., Aschner, J.L., Davis, P.G., McGrath-Morrow, S.A., Soll, R.F., et al., 2019. Bronchopulmonary dysplasia. Nat. Rev. Dis. Primers 5, 78. https://doi.org/10.1038/s41572-019-0127-7.

Türei, D., Korcsmáros, T., Saez-Rodriguez, J., 2016. OmniPath: guidelines and gateway for literature-curated signaling pathway resources. Nat. Methods 13, 966–967. https://doi.org/10.1038/nmeth.4077.

Vila Ellis, L., Bywaters, J.D., Ceas, A., Liu, Y., Sucre, J.M., Chen, J., 2025. P53 maintains lineage fidelity during lung capillary injury-repair in neonatal hyperoxia. JCI Insight 10, e182880. https://doi.org/10.1172/jci.insight.182880.

Wei, Z., Wang, Y., Gao, Y., Wang, S., Li, P., Si, D., Gao, Y., Wu, S., Li, D., Dong, K., et al., 2026. Benchmarking algorithms for generalizable single-cell perturbation response prediction. Nat. Methods 23, 451–464. https://doi.org/10.1038/s41592-025-02980-0.

Zanini, F., Che, X., Knutsen, C., Liu, M., Suresh, N.E., Domingo-Gonzalez, R., Dou, S.H., Zhang, D., Pryhuber, G.S., Jones, R.C., et al., 2023. Developmental diversity and unique sensitivity to injury of lung endothelial subtypes during postnatal growth. iScience 26, 106097. https://doi.org/10.1016/j.isci.2023.106097.

## Tables

### Table 1. Cross-cohort endothelial programs prioritized by the integrated analysis

| Program | Representative molecular evidence | Cross-cohort support | Biological interpretation |
|---|---|---|---|
| Mitochondrial bioenergetics | Oxidative phosphorylation, cellular respiration, ATP metabolism and mitochondrial translation | Directional concordance and FDR < 0.05 in all three designated external mouse contrasts | Sustained respiratory and energy-metabolic adaptation is a core endothelial response to neonatal hyperoxia |
| Oxidative stress and p53-associated injury | *Cdkn1a*, *Gdf15*, *Phlda3*, *Zmat3*, *Bax* and *Rps27l*; reactive-oxygen-species and NFE2L2 pathways | Listed genes belong to the 25-gene external meta-replicated set; ROS and NFE2L2 pathways passed FDR < 0.05 in all three external contrasts | A distributed injury network combines p53-associated transcription with antioxidant signaling |
| Matrix and wound remodeling | *Ecm1*, *Hapln1*, *Tinagl1*, *Serpine2* and *Lgals1*; EMT-like and wound-healing pathways | Gene directions reproduced across all external contrasts; both pathways passed FDR < 0.05 in all three contrasts | Structural remodeling accompanies metabolic and oxidative stress responses |
| Capillary and vascular signaling | *Apln*, *Ddah1* and *Cyp26b1*; Apln–Aplnr, Gdf15–Tgfbr2 and Tgfb2–Acvr1/Tgfbr1 routes | Gene directions reproduced across all external contrasts; ligand–receptor routes recurred as expression-supported hypotheses | Candidate signaling axes connect capillary maintenance, injury signaling and TGF-β-responsive remodeling |


## Figure legends

### Fig. 1. Study design links cross-cohort mouse biology to prediction and human concordance.

A: Dataset roles and analysis locks. Discovery and internal animal-level replication fixed gene and pathway sets before three external mouse cohorts, two independently trained prediction cohorts and the human atlas were evaluated. B: Each result class supports a distinct claim: reproducible mouse biology, cohort-sensitive prediction and a bounded human-transfer conclusion. EC, endothelial cells; BPD, bronchopulmonary dysplasia.

### Fig. 2. Neonatal hyperoxia induces a multi-subtype endothelial stress and remodeling program.

A: Purified P14 endothelial discovery effects. Dashed vertical lines mark absolute log2 fold-change of 1; outlined points are members of the 33-gene multi-subtype set. B: Discovery and internal animal-by-subtype log2 fold-changes for capillary and aerocyte-capillary endothelium; outlined points denote locked genes. C: Direction concordance among significant overlaps across five endothelial subtypes. Labels give concordant/evaluable genes. D: Selected directionally replicated pathway programs across discovery and internal endothelial analyses. Cap, general capillary; Cap-a, aerocyte capillary; Art, arterial; Lymph, lymphatic; NES, normalized enrichment score.

### Fig. 3. Independent mouse cohorts retain injury, bioenergetic and remodeling signals.

A: Discovery-oriented log2 fold-changes for 29 multi-subtype genes evaluable in all three external studies. B: Signed external meta-analysis effect versus false-discovery evidence; purple points are 33-gene-set members. C: Numbers of evaluable, directionally concordant and signed-meta-FDR-replicated genes for the 198- and 33-gene sets. D: Fractions of locked pathways directionally concordant in all three designated study-level contrasts. Full cohort-level effects are shown in Fig. S1. GO_BP, Gene Ontology Biological Process.

### Fig. 4. Perturbation-prediction rankings change across independently trained cohorts.

A: Held-out-age Spearman correlations for six GSE151974 capillary tasks using fold-specific 1,800-gene HVG endpoints. B: Paired joint-animal-bootstrap differences relative to PCA shift in GSE151974; points are median differences and bars are 95% percentile intervals. C: Corresponding four-task point estimates in the independent eight-animal GSE243129 cohort. D: Mean point estimates in the two cohorts. GSE243129 contains two animals per age-by-oxygen group; its ordering is hypothesis-generating. Exact animal-subset analyses are shown in Fig. S2. Cap, general capillary; Cap-a, aerocyte capillary; Art, arterial; Lymph, lymphatic; OT, optimal transport.

### Fig. 5. Human BPD-associated endothelium shows little concordance with the locked mouse program.

A: Donor-level adjusted coefficients and HC3 95% confidence intervals for strict 25-gene and extended 33-gene scores. Shading marks ±0.5 standard deviations. Gestational age and BPD status were completely confounded, so the coefficient represents the combined prematurity/BPD state. B: Mouse discovery effects versus adjusted human gene-level effects for the strict set. C: Full-pseudobulk, equal-cell and equal-UMI sensitivities. D: Fractions of 156 testable pathways matching the locked mouse direction; none met the joint direction-plus-human-FDR criterion. CAP1 and CAP2 denote the two deposited capillary populations.
