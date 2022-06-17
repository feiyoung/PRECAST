

Melt <- Seurat:::Melt
MultiExIPlot <- Seurat:::MultiExIPlot
SingleExIPlot <- Seurat:::SingleExIPlot
InvertHex <- Seurat:::InvertHex
"%||%" <- Seurat:::"%||%"

utils::globalVariables(c("Human_HK_genes", "Mouse_HK_genes",
        "makeCluster", "stopCluster", 
        "Melt", "MultiExIPlot", "SingleExIPlot", "InvertHex", "%||%"))
