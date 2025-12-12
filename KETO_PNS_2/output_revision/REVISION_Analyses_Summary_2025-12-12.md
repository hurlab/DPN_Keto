# Revision Analyses Summary (Dec 12, 2025)

## Overview
- Reran revision analyses in `output_revision` with updated outputs:
  - Figures: PDFs/PNGs in `output_revision/Figures/` (Top10/Top20 GO/KEGG heatmaps with and without values; Venns with larger labels).
  - Tables: `output_revision/Tables/` (DEGs/unique DEGs now include `Gene`, `log2FoldChange`, `pvalue`; overlap files also include per-comparison log2FC/pvalues).
- GO/KEGG use richR (live annotations, `builtin=FALSE`); all cutoffs unchanged from prior runs.

## Analysis 1 – DR-subtracted intervention effects (vs HFD)
- Method: Removed DRvsHFD DEGs from each intervention vs HFD; Venn + unique sets; GO/KEGG on uniques.
- Unique DEGs: EX 283; KDI 158; KDI_EX 285.
- KEGG highlights (unique sets):
  - EX: Inositol phosphate metabolism; O-glycan biosynthesis; Axon guidance.
  - KDI: ECM-receptor interaction; Focal adhesion; Cell adhesion molecules.
  - KDI_EX: Fatty acid degradation; Butanoate/propanoate metabolism; Peroxisome.
- Key files: `Analysis1_*` in Tables/Figures (e.g., `Analysis1_EXvsHFD_Unique_DEGs.csv`, `Analysis1_DR_Subtracted_Interventions_Venn.pdf/png`).

## Analysis 2 – Interventions vs SD, filtered against DR (dietary reversal) effects
- Method: From each intervention vs SD DEG set, removed all DRvsSD DEGs to isolate intervention-specific benefits beyond diet reversal; Venn + uniques; GO/KEGG on uniques.
- Unique DEGs (DR-filtered): EX 436; KDI 563; KDI_EX 621.
- KEGG highlights (DR-filtered uniques):
  - EXvsSD: Small cell lung cancer; Protein digestion/absorption; ECM/PI3K-Akt signaling (collagen-driven).
  - KDIvsSD: Hematopoietic cell lineage; Influenza A; Antigen presentation/immune signaling.
  - KDI_EXvsSD: Thyroid hormone signaling; mTOR signaling; AMPK/PI3K-Akt axis.
- Key files: prefixed `Analysis2_DRFiltered_*` (e.g., `Analysis2_DRFiltered_Interventions_Venn.pdf/png`, `Analysis2_EXvsSD_DRFiltered_Unique_DEGs.csv`, GO/KEGG CSVs).

## Analysis 3 – Nerve-specific intervention and maintenance schemes
- Method: For each comparison, removed matching gastroc DEGs from sciatic DEGs; Venns for intervention (EX/DR/KDI/KDI_EX/KD vs HFD) and maintenance (HFDvsSD, KDvsSD, KDvsHFD); GO/KEGG on nerve-specific uniques.
- Nerve-specific unique DEGs:
  - Intervention: DR 285; EX 283; KDI 167; KDI_EX 282; KD 197.
  - Maintenance: HFDvsSD 950; KDvsSD 123.
- KEGG highlights (examples):
  - EX nerve-specific: Inositol phosphate metabolism; O-glycan biosynthesis.
  - KDI nerve-specific: ECM-receptor interaction; Cell adhesion molecules.
  - KDI_EX nerve-specific: Fatty acid degradation; Butanoate/propanoate metabolism.
  - HFDvsSD nerve-specific: Tight junction; PI3K-Akt; Hippo signaling.
- Key files: `Analysis3_*` in Tables/Figures (e.g., `Analysis3_NerveSpecific_Interventions_Venn.pdf/png`, `Analysis3_EXvsHFD_NerveSpecific_Unique_DEGs.csv`, associated GO/KEGG CSVs).

## Notes on outputs
- DEG/Unique/overlap tables now include log2FC and p-value columns per comparison (padj removed to match DEG calling criteria).
- Figures exist in both PDF and PNG; heatmaps labeled with `-log10(<metric>)` and optional overlaid values.
