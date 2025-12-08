# Quick test of VennDiagram
library(VennDiagram)

# Create test gene sets
set1 <- paste0("Gene", 1:100)
set2 <- paste0("Gene", 50:150)
set3 <- paste0("Gene", 100:200)

# Create 3-way Venn
pdf("test_venn.pdf")
venn.plot <- venn.diagram(
  x = list(A = set1, B = set2, C = set3),
  category.names = c("Set A", "Set B", "Set C"),
  fill = c("red", "blue", "green")
)
grid.draw(venn.plot)
dev.off()

cat("Venn diagram test completed!\n")