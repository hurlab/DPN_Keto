# Revision Analyses Summary (Dec 12, 2025)

## Overview
- Reran revision analyses in `output_revision` with updated outputs:
  - Figures: PDFs/PNGs in `output_revision/Figures/` (Top10/Top20 GO/KEGG heatmaps with and without values; Venns with larger labels).
  - Tables: `output_revision/Tables/` (DEGs/unique DEGs now include `Gene`, `log2FoldChange`, `pvalue`; overlap files also include per-comparison log2FC/pvalues).
- GO/KEGG use richR (live annotations, `builtin=FALSE`); all cutoffs unchanged from prior runs.

## Analysis 1 – DR-subtracted intervention effects (vs HFD)
- Method: Removed DRvsHFD DEGs from each intervention vs HFD; Venn + unique sets; GO/KEGG on uniques.
- Unique DEGs: EX 283; KDI 158; KDI_EX 285.
- Interpretation:
  - EX (on HFD) shows neuronal/synaptic remodeling (axon guidance) and phosphoinositide signaling, indicating activity-driven nerve adaptations beyond diet withdrawal.
  - KDI enriches cell-adhesion/ECM signaling (ECM-receptor, CAMs), suggesting structural remodeling distinct from DR effects.
  - KDI_EX engages lipid catabolism (fatty acid degradation, peroxisome) and short-chain fatty acid pathways, pointing to combined dietary/exercise metabolic reprogramming.
- Key files: `Analysis1_*` in Tables/Figures (e.g., `Analysis1_EXvsHFD_Unique_DEGs.csv`, `Analysis1_DR_Subtracted_Interventions_Venn.pdf/png`).

## Analysis 2 – Interventions vs DR baseline
- Method: Direct contrasts EXvsDR, KDIvsDR, KDI_EXvsDR (DESeq2 from `scndds`); DEGs called at p < 0.05; Venn + uniques; GO/KEGG on uniques.
- Unique DEGs: EXvsDR 1,681; KDIvsDR 836; KDI_EXvsDR 832.
- Interpretation:
  - EXvsDR: strong cell-cycle/microtubule/motor protein programs (Cell cycle, Motor proteins), consistent with exercise-driven proliferative/repair signatures beyond DR.
  - KDIvsDR: translational and neurodegenerative modules (Ribosome, ALS) suggesting KD-specific proteostasis/mitochondrial adaptations not seen with DR.
  - KDI_EXvsDR: circadian and PI3K-Akt signaling (circadian rhythm, PI3K-Akt), indicating combined diet+exercise reprograms timing and growth-factor pathways beyond DR.
- Key files: prefixed `Analysis2_*` (e.g., `Analysis2_DR_Interventions_Venn.pdf/png`, `Analysis2_EXvsDR_Unique_DEGs.csv`, GO/KEGG CSVs, Top10/Top20 heatmaps under `Analysis2_DR_GO/KEGG_*`).

## Response to Reviewer 1 (Comment 2) – Intervention effects vs diet reversal
- Concern: Intervention benefits might reflect generic HFD withdrawal; need DR comparisons.
- Primary evidence (Analysis 1, DR-subtracted vs HFD): EX 283, KDI 158, KDI_EX 285 unique DEGs after removing DRvsHFD genes, showing synaptic/phosphoinositide (EX), ECM/CAM (KDI), and lipid catabolism/SCFA (KDI_EX) pathways—intervention-specific signals not explained by diet reversal.
- Auxiliary evidence (Analysis 2, interventions vs DR): direct DR baseline contrasts further show distinct intervention signatures (cell cycle/motor proteins for EX; ribosome/proteostasis for KDI; circadian/PI3K-Akt for KDI_EX).
- Deliverables: `Analysis1_*` (core DR-subtracted evidence) and `Analysis2_*` (auxiliary DR baseline) tables/figures with log2FC/pvalues.

## Analysis 3 – Nerve-specific intervention and maintenance schemes
- Method: For each comparison, removed matching gastroc DEGs from sciatic DEGs; Venns for intervention (EX/DR/KDI/KDI_EX/KD vs HFD) and maintenance (HFDvsSD, KDvsSD, KDvsHFD); GO/KEGG on nerve-specific uniques.
- Nerve-specific unique DEGs:
  - Intervention: DR 285; EX 283; KDI 167; KDI_EX 282; KD 197.
  - Maintenance: HFDvsSD 950; KDvsSD 123.
- Interpretation:
  - EX nerve-specific: phosphoinositide and glycan remodeling imply neuronal signaling/plasticity changes unique to exercise under HFD.
  - KDI nerve-specific: ECM/CAM enrichment suggests nerve structural remodeling distinct from muscle.
  - KDI_EX nerve-specific: fatty acid and short-chain fatty acid metabolism point to enhanced neural bioenergetics from combined intervention.
  - HFDvsSD nerve-specific: tight junction/PI3K-Akt/Hippo highlight barrier and growth-control shifts specific to nerve under metabolic stress.
- Key files: `Analysis3_*` in Tables/Figures (e.g., `Analysis3_NerveSpecific_Interventions_Venn.pdf/png`, `Analysis3_EXvsHFD_NerveSpecific_Unique_DEGs.csv`, associated GO/KEGG CSVs).

## Response to Reviewer 3 (nerve specificity)
- Concern: Need to demonstrate nerve (SCN) effects are not generic muscle responses.
- Action: Analysis 3 removed gastroc DEGs to retain nerve-only signals; Venns and enrichments computed on nerve-specific uniques.
- Evidence: Large nerve-only sets (e.g., HFDvsSD 950; EXvsHFD 283; DRvsHFD 285; KDI_EXvsHFD 282; KDIvsHFD 167; KDvsHFD 197; KDvsSD 123) with nerve-focused pathways (tight junction/PI3K-Akt/Hippo; phosphoinositide/axon guidance; peroxisome/PPAR; fatty acid degradation/SCFA; ECM/CAM; TCA/OXPHOS).
- Deliverables: `Analysis3_*` Venns/GO/KEGG files with log2FC/pvalues to document nerve-specific biology.

## Notes on outputs
- DEG/Unique/overlap tables now include log2FC and p-value columns per comparison (padj removed to match DEG calling criteria).
- Figures exist in both PDF and PNG; heatmaps labeled with `-log10(<metric>)` and optional overlaid values.
