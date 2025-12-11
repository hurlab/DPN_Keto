# Group name convention
# Note that there are some discrepancies in group naming between the manuscript, scripts, and review reports
# Use the following notion
# SDR = DR; KDI-EXI = KDI_EX; EXI = EX
 
Reviewer 1, Comment 2: From my understanding, a central interpretive gap is the absence of comparisons between intervention groups (KDI, EXI, KDI-EXI) and the SDI diet-reversal group. As currently presented, all intervention effects are interpreted only relative to HFD? This leaves open the possibility that some improvements are driven by stopping HFD rather than by the interventions themselves. Analyses of KDI vs SDI, EXI vs SDI, and KDI-EXI vs SDI could help determine intervention-specific benefits. Authors could at least perform an analysis that removes DEGs shared between SDI vs HFD and each intervention vs HFD. This would isolate genes and pathways uniquely affected by KDI, EXI, or KDI-EXI rather than those reflecting generic reversal of HFD-induced injury.

=> Analysis #1: There is a four-way Venn Diagram analysis performed with four DEG sets Figure 3B. 
   DR vs HFD
   EX vs HFD
   KDI vs HFD
   KDI_EX vs HFD
   
   First remove all DR-vs-HFD DEGs from the remaining three sets, and then perform 3-way Venn diagram using VennDetail package. This will give the DEGs not related to DR vs HFD (stopping HFD). Then, for each unique subset from this 3-way VennDiagram, meaning only those DEGs belonging to the comparison and not shared by other comparisons. There should be three such DEG sets, which should be 283, 285, and 158 genes. Perform functional enrichment anlayses using richR package with GO and KEGG. Create summary heatmaps (one with GO and another with KEGG) using, top 10 or top 20 most significantly enrichend annotations, visualizing -log10(pvalue) as the color. these heatmap will collectively show the top functinos as well as their overlap visually. 

==> Analysis #2: As suggested by the reviewer, let's look at "KDI vs SDI, EXI vs SDI, and KDI-EXI vs SDI" this three-way comparison. Same technical approach as above, including Venn diagram, GO/KEGG enrichment, heatmap. 


Reviewer 3, Comment 7: I suspect that many of the gene changes are related to correction of metabolic syndrome by treatments with ketogenic diet or exercise.  Could the muscle be used as a control for the changes in the nerve. Authors should compare gene changes in muscle to those in nerve to assess for this. Would be particularly interested in gene changes that are happening in nerve and NOT in muscle.

==> Analysis #3: This basically ask for comparisons between nerve (sciatic) vs muscle (gastroc). In addition to the 'intervention' scheme mentioned above with 4 DEG sets, there is 'maintenance' scheme with three DEG sets, including HFD vs SD, KD vs SD, and KD vs HFD. For each comparison (like HFD vs SD), let's substract gastroc DEGs from sciatic nerve DEGs, so that we can keep only the nerve DEGs. Then, we can do additional venn diagram analyses for 'intervention' - nerve-specific and 'maintenance' - nerve-specific. Functional enrichments and heatmaps can also be generated . 


