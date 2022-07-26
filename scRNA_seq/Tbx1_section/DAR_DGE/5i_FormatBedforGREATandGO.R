### Taking motif filtered DARs per cluster and formating for input into GREAT and GO
### First GREAT:
setwd("/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/motif_peaks")

## AHF
AHF_motif1 <- read.csv(file="AHF/AHF_motif1_outputfile.txt", sep = '\t')
AHF_motif2 <- read.csv(file="AHF/AHF_motif2_outputfile.txt", sep = '\t')
AHF_motif3 <- read.csv(file="AHF/AHF_motif3_outputfile.txt", sep = '\t')
AHF_motif4 <- read.csv(file="AHF/AHF_motif4_outputfile.txt", sep = '\t')
AHF_motif5 <- read.csv(file="AHF/AHF_motif5_outputfile.txt", sep = '\t')
AHF_motif6 <- read.csv(file="AHF/AHF_motif6_outputfile.txt", sep = '\t')

# filter rows without a motif
AHF_motif1_filtered <- AHF_motif1[AHF_motif1$X1.GTGWTGACAGAY.BestGuess.Tbx20.T.box..Heart.Tbx20.ChIP.Seq.GSE29636..Homer.0.860..Distance.From.Peak.sequence.strand.conservation.!="",]
AHF_motif2_filtered <- AHF_motif2[AHF_motif2$X7.GACARACG.BestGuess.ATOH7.MA1468.1.Jaspar.0.716..Distance.From.Peak.sequence.strand.conservation. !="",]
AHF_motif3_filtered <- AHF_motif3[AHF_motif3$X8.CAGCTGTTTK.BestGuess.BHLHA15.var.2..MA1472.1.Jaspar.0.947..Distance.From.Peak.sequence.strand.conservation.!="",]
AHF_motif4_filtered <- AHF_motif4[AHF_motif4$X8.BYTAATTR.BestGuess.PDX1.MA0132.2.Jaspar.0.951..Distance.From.Peak.sequence.strand.conservation.!="",]
AHF_motif5_filtered <- AHF_motif5[AHF_motif5$X10.STBTGTTTAY.BestGuess.FOXL1.MA0033.2.Jaspar.0.904..Distance.From.Peak.sequence.strand.conservation.!="",]
AHF_motif6_filtered <- AHF_motif6[AHF_motif6$X11.CCTTATCT.BestGuess.GATA4.MA0482.2.Jaspar.0.987..Distance.From.Peak.sequence.strand.conservation.!="",]



## Format for GREAT
AHF_motif1_filtered_great <- AHF_motif1_filtered[,c(2,3,4,1,6,5)]
AHF_motif2_filtered_great <- AHF_motif2_filtered[,c(2,3,4,1,6,5)]
AHF_motif3_filtered_great <- AHF_motif3_filtered[,c(2,3,4,1,6,5)]
AHF_motif4_filtered_great <- AHF_motif4_filtered[,c(2,3,4,1,6,5)]
AHF_motif5_filtered_great <- AHF_motif5_filtered[,c(2,3,4,1,6,5)]
AHF_motif6_filtered_great <- AHF_motif6_filtered[,c(2,3,4,1,6,5)]

# save files
write.table(AHF_motif1_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/AHF/AHF_motif1_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')
write.table(AHF_motif2_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/AHF/AHF_motif2_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')
write.table(AHF_motif3_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/AHF/AHF_motif3_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')
write.table(AHF_motif4_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/AHF/AHF_motif4_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')
write.table(AHF_motif5_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/AHF/AHF_motif5_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')
write.table(AHF_motif6_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/AHF/AHF_motif6_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')


#### 
## Pharyngeal Mesoderm
PharyngealMesoderm_motif1 <- read.csv(file="PharyngealMesoderm/PharyngealMesoderm_motif1_outputfile.txt", sep = '\t')
PharyngealMesoderm_motif2 <- read.csv(file="PharyngealMesoderm/PharyngealMesoderm_motif2_outputfile.txt", sep = '\t')
PharyngealMesoderm_motif3 <- read.csv(file="PharyngealMesoderm/PharyngealMesoderm_motif3_outputfile.txt", sep = '\t')
PharyngealMesoderm_motif4 <- read.csv(file="PharyngealMesoderm/PharyngealMesoderm_motif4_outputfile.txt", sep = '\t')
PharyngealMesoderm_motif5 <- read.csv(file="PharyngealMesoderm/PharyngealMesoderm_motif5_outputfile.txt", sep = '\t')
PharyngealMesoderm_motif6 <- read.csv(file="PharyngealMesoderm/PharyngealMesoderm_motif6_outputfile.txt", sep = '\t')

# filter rows without a motif
PharyngealMesoderm_motif1_filtered <- PharyngealMesoderm_motif1[PharyngealMesoderm_motif1$X1.GTGWTGACAGAT.BestGuess.Tbx20.T.box..Heart.Tbx20.ChIP.Seq.GSE29636..Homer.0.885..Distance.From.Peak.sequence.strand.conservation. !="",]
PharyngealMesoderm_motif2_filtered <- PharyngealMesoderm_motif2[PharyngealMesoderm_motif2$X6.CTGTTAAY.BestGuess.PB0109.1_Bbx_2.Jaspar.0.822..Distance.From.Peak.sequence.strand.conservation. !="",]
PharyngealMesoderm_motif3_filtered <- PharyngealMesoderm_motif3[PharyngealMesoderm_motif3$X5.CAGCTGBYNN.BestGuess.Twist2.bHLH..Myoblast.Twist2.Ty1.ChIP.Seq.GSE127998..Homer.0.964..Distance.From.Peak.sequence.strand.conservation. !="",]
PharyngealMesoderm_motif4_filtered <- PharyngealMesoderm_motif4[PharyngealMesoderm_motif4$X6.BTAATTAGMN.BestGuess.DLX2.Homeobox..BasalGanglia.Dlx2.ChIP.seq.GSE124936..Homer.0.977..Distance.From.Peak.sequence.strand.conservation. !="",]
PharyngealMesoderm_motif5_filtered <- PharyngealMesoderm_motif5[PharyngealMesoderm_motif5$X6.TCCCCAGGGRVN.BestGuess.EBF1.EBF..Near.E2A.ChIP.Seq.GSE21512..Homer.0.957..Distance.From.Peak.sequence.strand.conservation. !="",]
PharyngealMesoderm_motif6_filtered <- PharyngealMesoderm_motif6[PharyngealMesoderm_motif6$X8.YWNKCTGHCA.BestGuess.Meis1.Homeobox..MastCells.Meis1.ChIP.Seq.GSE48085..Homer.0.699..Distance.From.Peak.sequence.strand.conservation. !="",]

## Format for GREAT
PharyngealMesoderm_motif1_filtered_great <- PharyngealMesoderm_motif1_filtered[,c(2,3,4,1,6,5)]
PharyngealMesoderm_motif2_filtered_great <- PharyngealMesoderm_motif2_filtered[,c(2,3,4,1,6,5)]
PharyngealMesoderm_motif3_filtered_great <- PharyngealMesoderm_motif3_filtered[,c(2,3,4,1,6,5)]
PharyngealMesoderm_motif4_filtered_great <- PharyngealMesoderm_motif4_filtered[,c(2,3,4,1,6,5)]
PharyngealMesoderm_motif5_filtered_great <- PharyngealMesoderm_motif5_filtered[,c(2,3,4,1,6,5)]
PharyngealMesoderm_motif6_filtered_great <- PharyngealMesoderm_motif6_filtered[,c(2,3,4,1,6,5)]

# save files
write.table(PharyngealMesoderm_motif1_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/PharyngealMesoderm/PharyngealMesoderm_motif1_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')
write.table(PharyngealMesoderm_motif2_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/PharyngealMesoderm/PharyngealMesoderm_motif2_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')
write.table(PharyngealMesoderm_motif3_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/PharyngealMesoderm/PharyngealMesoderm_motif3_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')
write.table(PharyngealMesoderm_motif4_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/PharyngealMesoderm/PharyngealMesoderm_motif4_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')
write.table(PharyngealMesoderm_motif5_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/PharyngealMesoderm/PharyngealMesoderm_motif5_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')
write.table(PharyngealMesoderm_motif6_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/PharyngealMesoderm/PharyngealMesoderm_motif6_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')

## Endoderm
Endoderm_motif1 <- read.csv(file="Endoderm/Endoderm_motif1_outputfile.txt", sep = '\t')
Endoderm_motif2 <- read.csv(file="Endoderm/Endoderm_motif2_outputfile.txt", sep = '\t')
Endoderm_motif3 <- read.csv(file="Endoderm/Endoderm_motif3_outputfile.txt", sep = '\t')
Endoderm_motif4 <- read.csv(file="Endoderm/Endoderm_motif4_outputfile.txt", sep = '\t')
Endoderm_motif5 <- read.csv(file="Endoderm/Endoderm_motif5_outputfile.txt", sep = '\t')
Endoderm_motif6 <- read.csv(file="Endoderm/Endoderm_motif6_outputfile.txt", sep = '\t')

# filter rows without a motif
Endoderm_motif1_filtered <- Endoderm_motif1[Endoderm_motif1$X1.GTGHTGACAGAY.BestGuess.Tbx20.T.box..Heart.Tbx20.ChIP.Seq.GSE29636..Homer.0.849..Distance.From.Peak.sequence.strand.conservation. !="",]
Endoderm_motif2_filtered <- Endoderm_motif2[Endoderm_motif2$X4.GAAGATTGATGG.BestGuess.HOXA1.Homeobox..mES.Hoxa1.ChIP.Seq.SRP084292..Homer.0.815..Distance.From.Peak.sequence.strand.conservation. !="",]
Endoderm_motif3_filtered <- Endoderm_motif3[Endoderm_motif3$X6.CAATTAGC.BestGuess.EN2.MA0642.1.Jaspar.0.913..Distance.From.Peak.sequence.strand.conservation. !="",]
Endoderm_motif4_filtered <- Endoderm_motif4[Endoderm_motif4$X7.CAAACAAAGGAG.BestGuess.Sox3.MA0514.1.Jaspar.0.876..Distance.From.Peak.sequence.strand.conservation. !="",]
Endoderm_motif5_filtered <- Endoderm_motif5[Endoderm_motif5$X8.TGTTAACACA.BestGuess.PB0109.1_Bbx_2.Jaspar.0.805..Distance.From.Peak.sequence.strand.conservation. !="",]
Endoderm_motif6_filtered <- Endoderm_motif6[Endoderm_motif6$X9.AAGTGCTTTG.BestGuess.Nkx3.1.Homeobox..LNCaP.Nkx3.1.ChIP.Seq.GSE28264..Homer.0.762..Distance.From.Peak.sequence.strand.conservation. !="",]

## Format for GREAT
Endoderm_motif1_filtered_great <- Endoderm_motif1_filtered[,c(2,3,4,1,6,5)]
Endoderm_motif2_filtered_great <- Endoderm_motif2_filtered[,c(2,3,4,1,6,5)]
Endoderm_motif3_filtered_great <- Endoderm_motif3_filtered[,c(2,3,4,1,6,5)]
Endoderm_motif4_filtered_great <- Endoderm_motif4_filtered[,c(2,3,4,1,6,5)]
Endoderm_motif5_filtered_great <- Endoderm_motif5_filtered[,c(2,3,4,1,6,5)]
Endoderm_motif6_filtered_great <- Endoderm_motif6_filtered[,c(2,3,4,1,6,5)]

# save files
write.table(Endoderm_motif1_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/Endoderm/Endoderm_motif1_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')
write.table(Endoderm_motif2_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/Endoderm/Endoderm_motif2_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')
write.table(Endoderm_motif3_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/Endoderm/Endoderm_motif3_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')
write.table(Endoderm_motif4_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/Endoderm/Endoderm_motif4_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')
write.table(Endoderm_motif5_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/Endoderm/Endoderm_motif5_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')
write.table(Endoderm_motif6_filtered_great, file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GREAT_beds/Endoderm/Endoderm_motif6_filtered_great.bed", row.names = F,col.names = F, quote = F, sep = '\t')


####################################
# Output for GO

## AHF
AHF_motif1_filtered_GO <- as.factor(AHF_motif1_filtered[,16])
AHF_motif2_filtered_GO <- as.factor(AHF_motif2_filtered[,16])
AHF_motif3_filtered_GO <- as.factor(AHF_motif3_filtered[,16])
AHF_motif4_filtered_GO <- as.factor(AHF_motif4_filtered[,16])
AHF_motif5_filtered_GO <- as.factor(AHF_motif5_filtered[,16])
AHF_motif6_filtered_GO <- as.factor(AHF_motif6_filtered[,16])

write.table(AHF_motif1_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/AHF/AHF_motif1_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')
write.table(AHF_motif2_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/AHF/AHF_motif2_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')
write.table(AHF_motif3_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/AHF/AHF_motif3_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')
write.table(AHF_motif4_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/AHF/AHF_motif4_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')
write.table(AHF_motif5_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/AHF/AHF_motif5_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')
write.table(AHF_motif6_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/AHF/AHF_motif6_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')


## PharyngealMesoderm
PharyngealMesoderm_motif1_filtered_GO <- as.factor(PharyngealMesoderm_motif1_filtered[,16])
PharyngealMesoderm_motif2_filtered_GO <- as.factor(PharyngealMesoderm_motif2_filtered[,16])
PharyngealMesoderm_motif3_filtered_GO <- as.factor(PharyngealMesoderm_motif3_filtered[,16])
PharyngealMesoderm_motif4_filtered_GO <- as.factor(PharyngealMesoderm_motif4_filtered[,16])
PharyngealMesoderm_motif5_filtered_GO <- as.factor(PharyngealMesoderm_motif5_filtered[,16])
PharyngealMesoderm_motif6_filtered_GO <- as.factor(PharyngealMesoderm_motif6_filtered[,16])

write.table(PharyngealMesoderm_motif1_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/PharyngealMesoderm/PharyngealMesoderm_motif1_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')
write.table(PharyngealMesoderm_motif2_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/PharyngealMesoderm/PharyngealMesoderm_motif2_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')
write.table(PharyngealMesoderm_motif3_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/PharyngealMesoderm/PharyngealMesoderm_motif3_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')
write.table(PharyngealMesoderm_motif4_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/PharyngealMesoderm/PharyngealMesoderm_motif4_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')
write.table(PharyngealMesoderm_motif5_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/PharyngealMesoderm/PharyngealMesoderm_motif5_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')
write.table(PharyngealMesoderm_motif6_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/PharyngealMesoderm/PharyngealMesoderm_motif6_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')

## Endoderm
Endoderm_motif1_filtered_GO <- as.factor(Endoderm_motif1_filtered[,16])
Endoderm_motif2_filtered_GO <- as.factor(Endoderm_motif2_filtered[,16])
Endoderm_motif3_filtered_GO <- as.factor(Endoderm_motif3_filtered[,16])
Endoderm_motif4_filtered_GO <- as.factor(Endoderm_motif4_filtered[,16])
Endoderm_motif5_filtered_GO <- as.factor(Endoderm_motif5_filtered[,16])
Endoderm_motif6_filtered_GO <- as.factor(Endoderm_motif6_filtered[,16])

write.table(Endoderm_motif1_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/Endoderm/Endoderm_motif1_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')
write.table(Endoderm_motif2_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/Endoderm/Endoderm_motif2_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')
write.table(Endoderm_motif3_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/Endoderm/Endoderm_motif3_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')
write.table(Endoderm_motif4_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/Endoderm/Endoderm_motif4_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')
write.table(Endoderm_motif5_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/Endoderm/Endoderm_motif5_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')
write.table(Endoderm_motif6_filtered_GO,file = "/Users/sranade/scATAC-seq/Tbx1_E925/DAR_homer/GO_csv/Endoderm/Endoderm_motif6_filtered_GO.csv", row.names = F, quote = F, col.names = F, sep = '\t')




