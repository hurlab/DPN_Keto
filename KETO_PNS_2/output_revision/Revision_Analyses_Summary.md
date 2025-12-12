# Revision Analyses Summary (PNS, R 4.5.2 bioinfo env)

## Overview
This run implements the three requested analyses and writes outputs to `output_revision/Figures` and `output_revision/Tables/Pathway_Analysis`. KEGG/GO enrichments use richR with live annotations (`builtin=FALSE`).

## Analysis 1 – DR-subtracted intervention effects (EX vs HFD, KDI vs HFD, KDI_EX vs HFD)
- **Method:** Removed all DRvsHFD DEGs from each intervention vs HFD set, then 3-way Venn (EX, KDI, KDI_EX).
- **Highlight Unique DEGs:** EX 283; KDI 158; KDI_EX 285
- **Top KEGG (unique sets):**
  - EX: Inositol phosphate metabolism; Axon guidance; Glutamatergic synapse.
  - KDI: CAM interaction; ECM-receptor interaction; Focal adhesion.
  - KDI_EX: Butanoate metabolism; Fatty acid degradation; Propanoate metabolism.
- **Top GO (unique sets):**
  - EX: Dendritic spine organization; Regulation of synapse structure or activity.
  - KDI: Regulation of synapse structure or activity; Cell junction organization.
  - KDI_EX: Positive regulation of smooth muscle cell proliferation; Muscle tissue development.
- **Outputs:** `Analysis1_*` CSVs and Venn/heatmaps in `Figures/`.

## Analysis 2 – Interventions vs SD (EX vs SD, KDI vs SD, KDI_EX vs SD)
- **Method:** 3-way Venn across SD baselines.
- **Unique DEGs:** EX 455; KDI 624; KDI_EX 683.
- **Top KEGG (unique sets):**
  - EX: Small cell lung cancer; PI3K-Akt signaling; Ras signaling.
  - KDI: Human T-cell leukemia virus 1 infection; Th17 differentiation; Chemokine signaling.
  - KDI_EX: Thyroid hormone signaling; AMPK signaling; Insulin resistance.
- **Outputs:** `Analysis2_*` CSVs and Venn/heatmaps in `Figures/`.

## Analysis 3 – Nerve-specific DEGs (sciatic minus gastroc) for intervention & maintenance schemes
- **Method:** For each comparison, removed gastroc DEGs from sciatic DEGs; separate Venns for intervention (DR/EX/KDI/KDI_EX vs HFD) and maintenance (HFDvsSD, KDvsSD, KDvsHFD).
- **Nerve-specific unique DEGs (selected):**
  - HFDvsSD: 950
  - EXvsHFD: 283; DRvsHFD: 285; KDIvsHFD: 167; KDI_EXvsHFD: 282; KDvsHFD: 197
  - KDvsSD: 123
- **Top KEGG (nerve-specific uniques):**
  - HFDvsSD: Hippo signaling; Tight junction; PI3K-Akt signaling.
  - EXvsHFD: Inositol phosphate metabolism; Axon guidance.
  - DRvsHFD: Peroxisome; PPAR signaling.
  - KDI_EXvsHFD: Butanoate metabolism; Fatty acid degradation.
  - KDIvsHFD: CAM interaction; ECM-receptor interaction.
  - KDvsHFD: Citrate cycle (TCA); Oxidative phosphorylation.
- **Outputs:** `Analysis3_*` CSVs and Venn/heatmaps in `Figures/`.

## Point-by-point responses to reviewer comments
### Reviewer 1, Comment 2 (intervention effects beyond HFD reversal)
- We performed the DR-subtracted comparison: after removing DRvsHFD DEGs, the Venn shows EX 283, KDI 158, and KDI_EX 285 unique genes, demonstrating intervention-specific regulation distinct from diet reversal. KEGG highlights show synaptic/axon guidance (EX), cell-adhesion/ECM (KDI), and short-chain fatty acid metabolism (KDI_EX). GO highlights include synapse organization (EX, KDI) and smooth/muscle development (KDI_EX).
- We also compared interventions directly to SD (EXvsSD, KDIvsSD, KDI_EXvsSD). Unique gene sets are EX 455, KDI 624, and KDI_EX 683, indicating non-overlapping signatures. KEGG terms include PI3K-Akt/Ras (EX), immune/chemokine/Th17 (KDI), and thyroid hormone/AMPK/insulin-resistance (KDI_EX), supporting distinct biological effects relative to SD.
- Supporting outputs: `output_revision/Figures/Analysis1_*` and `Analysis2_*`; KEGG/GO tables in `output_revision/Tables/Pathway_Analysis/Analysis1_*` and `Analysis2_*`.

### Reviewer 3, Comment 7 (nerve specificity vs muscle)
- We subtracted gastroc DEGs from sciatic DEGs for both intervention and maintenance schemes. Nerve-specific unique DEGs remain substantial (e.g., HFDvsSD 950; EXvsHFD 283; DRvsHFD 285; KDIvsHFD 167; KDI_EXvsHFD 282; KDvsHFD 197; KDvsSD 123), confirming that sciatic nerve responses are largely not shared with muscle.
- KEGG enrichment of the nerve-only sets highlights neuronal and metabolic pathways: Hippo/tight junction (HFDvsSD); inositol phosphate/axon guidance (EXvsHFD); peroxisome/PPAR (DRvsHFD); butanoate/fatty acid degradation (KDI_EXvsHFD); CAM/ECM (KDIvsHFD); and TCA/OXPHOS (KDvsHFD), demonstrating nerve-specific biology rather than generalized muscle effects.
- Supporting outputs: `output_revision/Figures/Analysis3_*` and `output_revision/Tables/Pathway_Analysis/Analysis3_*`.

## Key Deliverables (paths relative to `KETO_PNS_2/output_revision`)
- Figures: `Figures/Analysis1_*_Venn.pdf`, `Analysis1_GO_Heatmap.pdf`, `Analysis1_KEGG_Heatmap.pdf`, analogous files for Analysis2 and Analysis3 (Intervention/Maintenance).
- Tables: `Tables/Pathway_Analysis/Analysis*_...csv` (filtered, unique DEGs, overlaps, KEGG/GO results).

## Notes
- GO enrichment now runs via `richGO` (no `builtin` flag). KEGG/GO annotations are pulled online (`builtin=FALSE`).
- All counts and top pathways come directly from the new CSVs produced in this run.
