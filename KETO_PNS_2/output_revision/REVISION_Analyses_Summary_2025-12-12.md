# Revision Analyses Summary (Dec 12, 2025)

## Overview
- Reran revision analyses in `output_revision` with updated outputs:
  - Figures: PDFs/PNGs in `output_revision/Figures/` (Top10/Top20 GO/KEGG heatmaps with and without values; Venns with larger labels).
  - Tables: `output_revision/Tables/` (DEGs/unique DEGs now include `Gene`, `log2FoldChange`, `pvalue`; overlap files also include per-comparison log2FC/pvalues).
- GO/KEGG use richR (live annotations, `builtin=FALSE`); all cutoffs unchanged from prior runs.

## Point-by-Point Responses (Reviewer 1, Comment 2)
- Concern: Intervention effects could be generic reversal of HFD; need comparisons against dietary reversal (DR/SD) to show intervention-specific benefits.
- Action: Added DR-filtered baseline analysis (Analysis 2) removing DRvsSD genes from EX/KDI/KDI_EX vs SD; this isolates intervention-specific DEGs beyond diet reversal.
- Evidence:
  - EX (still on HFD) retains 436 DR-filtered unique DEGs vs SD; top pathways (small cell lung cancer, ECM/PI3K-Akt) indicate extracellular remodeling/signaling not explained by stopping HFD.
  - KDI and KDI_EX show immune/endocrine programs (hematopoietic lineage, antigen presentation; thyroid hormone and mTOR signaling) absent from DR-only reversal.
  - DR-subtracted HFD contrasts (Analysis 1) corroborate intervention-unique genes vs HFD; EX retains synaptic/metabolic pathways despite HFD, reinforcing non-dietary effects.
- Deliverables: DR-filtered Venns/heatmaps (`Analysis2_DRFiltered_*`), and tables with log2FC/pvalues for all overlaps/uniques.

## Analysis 1 – DR-subtracted intervention effects (vs HFD)
- Method: Removed DRvsHFD DEGs from each intervention vs HFD; Venn + unique sets; GO/KEGG on uniques.
- Unique DEGs: EX 283; KDI 158; KDI_EX 285.
- Interpretation:
  - EX (on HFD) shows neuronal/synaptic remodeling (axon guidance) and phosphoinositide signaling, indicating activity-driven nerve adaptations beyond diet withdrawal.
  - KDI enriches cell-adhesion/ECM signaling (ECM-receptor, CAMs), suggesting structural remodeling distinct from DR effects.
  - KDI_EX engages lipid catabolism (fatty acid degradation, peroxisome) and short-chain fatty acid pathways, pointing to combined dietary/exercise metabolic reprogramming.
- Key files: `Analysis1_*` in Tables/Figures (e.g., `Analysis1_EXvsHFD_Unique_DEGs.csv`, `Analysis1_DR_Subtracted_Interventions_Venn.pdf/png`).

## Analysis 2 – Interventions vs SD, filtered against DR (dietary reversal) effects
- Method: From each intervention vs SD DEG set, removed all DRvsSD DEGs to isolate intervention-specific benefits beyond diet reversal; Venn + uniques; GO/KEGG on uniques.
- Unique DEGs (DR-filtered): EX 436; KDI 563; KDI_EX 621.
- Interpretation:
  - EXvsSD (still on HFD) retains ECM/PI3K-Akt signatures (collagen-rich “protein digestion/absorption”) implying EX-specific matrix/signaling shifts not achieved by DR.
  - KDIvsSD shows strong immune/antigen-presentation (hematopoietic lineage, antigen processing), indicating KD-driven immune modulation beyond diet reversal.
  - KDI_EXvsSD activates endocrine/metabolic control (thyroid hormone, mTOR, AMPK/PI3K-Akt), consistent with synergistic dietary+exercise benefits on energy signaling.
- Key files: prefixed `Analysis2_DRFiltered_*` (e.g., `Analysis2_DRFiltered_Interventions_Venn.pdf/png`, `Analysis2_EXvsSD_DRFiltered_Unique_DEGs.csv`, GO/KEGG CSVs).

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

## Notes on outputs
- DEG/Unique/overlap tables now include log2FC and p-value columns per comparison (padj removed to match DEG calling criteria).
- Figures exist in both PDF and PNG; heatmaps labeled with `-log10(<metric>)` and optional overlaid values.
