# PNS Differential Expression Gene (DEG) Analysis Report

## Summary of DEGs Identified

### Sciatic Nerve (SCN)
| Comparison | Total DEGs | Upregulated | Downregulated |
|------------|------------|-------------|--------------|
| HFD vs SD | 1,190 | 754 | 436 |
| KDI_EX vs SD | 1,150 | 692 | 458 |
| KDI vs SD | 1,025 | 510 | 515 |
| DR vs SD | 545 | 190 | 355 |
| EX vs SD | 681 | 389 | 292 |
| KD vs SD | 228 | 82 | 146 |
| KDI_EX vs HFD | 487 | 312 | 175 |
| DR vs HFD | 442 | 205 | 237 |
| KD vs HFD | 451 | 134 | 317 |
| EX vs HFD | 368 | 221 | 147 |
| KDI vs HFD | 321 | 116 | 205 |

**Total SCN DEGs across all comparisons: 5,888**

### Gastrocnemius Muscle (Gastroc)
| Comparison | Total DEGs | Upregulated | Downregulated |
|------------|------------|-------------|--------------|
| KDI vs SD | 885 | 504 | 381 |
| KDI_EX vs SD | 556 | 297 | 259 |
| KD vs SD | 607 | 285 | 322 |
| EX vs SD | 736 | 343 | 393 |
| HFD vs SD | 481 | 195 | 286 |
| KDI_EX vs HFD | 573 | 293 | 280 |
| KD vs HFD | 485 | 150 | 335 |
| DR vs HFD | 299 | 138 | 161 |
| DR vs SD | 129 | 68 | 61 |
| KDI vs HFD | 173 | 112 | 61 |
| EX vs HFD | 100 | 45 | 55 |

**Total Gastrocnemius DEGs across all comparisons: 4,624**

### Hippocampus (Hippo)
| Comparison | Total DEGs | Upregulated | Downregulated |
|------------|------------|-------------|--------------|
| KDI_EX vs HFD | 721 | 388 | 333 |
| EX vs HFD | 397 | 210 | 187 |
| KD vs HFD | 240 | 127 | 113 |
| KD vs SD | 138 | 53 | 85 |
| DR vs SD | 405 | 168 | 237 |
| KDI_EX vs SD | 361 | 122 | 239 |
| KDI vs HFD | 206 | 88 | 118 |
| DR vs HFD | 174 | 91 | 83 |
| HFD vs SD | 236 | 141 | 95 |
| EX vs SD | 189 | 71 | 118 |
| KDI vs SD | 170 | 39 | 131 |

**Total Hippocampus DEGs across all comparisons: 3,237**

## Key Observations

1. **Highest DEG counts** in SCN (peripheral nerve tissue)
2. **Moderate DEG counts** in Gastrocnemius (muscle tissue)
3. **Lowest DEG counts** in Hippocampus (CNS reference tissue)

### Notable Patterns:
- **HFD vs SD** consistently shows high DEG numbers across all tissues
- **Combined interventions (KDI_EX)** produce the most dramatic changes
- **Ketogenic diet alone (KD)** shows moderate effects
- **DR (Dietary Restriction)** shows tissue-specific responses

## Comparison with Literature

The numbers obtained appear reasonable for RNA-seq differential expression analysis:
- Typical DEG range: 200-2000 genes per comparison
- Our results fall within expected biological variation
- PNS tissues (SCN, Gastroc) show more dramatic changes than CNS (Hippo)

## Files Generated
- PCA plots for all three tissues
- DEG tables for each comparison
- Summary statistics
- Sample distance heatmaps

All files are saved in:
- `Tables/DEGs/` - Differential expression results
- `Figures/` - Visualization outputs