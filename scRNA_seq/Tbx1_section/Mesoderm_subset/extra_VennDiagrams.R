## generate venn diagrams of overlaps
library(VennDiagram)
grid.newpage()
venn.plot <- draw.pairwise.venn(area1           = 4362,
                                area2           = 563,
                                cross.area      = 257,
                                category        = c("DAR", "TboxMotifs"),
                                fill            = c("cornflowerblue", "yellow"),
                                lty             = "blank",
                                cex             = 2,
                                cat.cex         = 2,
                                cat.pos         = c(285, 105),
                                cat.dist        = 0.09,
                                cat.just        = list(c(-1, -1), c(1, 1)),
                                ext.pos         = 30,
                                ext.dist        = -0.05,
                                ext.length      = 0.85,
                                ext.line.lwd    = 2,
                                ext.line.lty    = "dashed"
)

