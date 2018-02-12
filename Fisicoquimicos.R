##### Script Ingrid Paola Figueroa #####

library(plotrix)
library(agricolae)

getwd()  # print the current working directory
setwd("~RamosFigueroa/Documents/OneDrive/Work/PROYECTO\ MANGLAR/BIOINF_ANALYSIS_MANGROVE_16S/Segundo\ Muestreo\ 16S/CanjerejosCamarones/Resultados\ Paola\ Figueroa/")



fisicoq <- read.delim("metadata.txt"); fisicoq <- as.data.frame(fisicoq)
fisicoq_promedios1 <- aggregate(fisicoq[, 5:25], list(fisicoq$groups), mean); row.names(fisicoq_promedios1) <- fisicoq_promedios1[,1]
fisicoq_promedios <- round(fisicoq_promedios1[, 2:22], 2)

fisicoq_errorstand1 <- aggregate(fisicoq[, 5:25], list(fisicoq$groups), std.error); row.names(fisicoq_errorstand1) <- fisicoq_errorstand1[,1]
fisicoq_errorstand <- round(fisicoq_errorstand1[, 2:22], 2)



write.table(fisicoq_promedios, file = "fisicoq_promedios.txt", sep = "\t", quote = FALSE, col.names = TRUE)
write.table(fisicoq_errorstand, file = "fisicoq_errorstand.txt", sep = "\t", quote = FALSE, col.names = TRUE)

#write.table(fisicoq_promediosT, file = "/Users/RamosFigueroa/Documents/OneDrive/Work/Tesis\ Maestr??a\ Microbiolog??a/Documentos\ de\ Tesis/fisicoq_promediosT.txt", sep = "\t", quote = FALSE, col.names = TRUE)
#write.table(fisicoq_errorstandT, file = "/Users/RamosFigueroa/Documents/OneDrive/Work/Tesis\ Maestr??a\ Microbiolog??a/Documentos\ de\ Tesis/fisicoq_errorstandT.txt", sep = "\t", quote = FALSE, col.names = TRUE)


# ANOVAS Y TUKEY

model_pH<-aov(pH~groups, data=fisicoq); summary(model_pH)
out_pH <- HSD.test(model_pH,"groups", group=TRUE,console=TRUE)
write.table(out_pH$groups, file = "out_pH.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_CE<-aov(CE~groups, data=fisicoq); summary(model_CE)
out_CE <- HSD.test(model_CE,"groups", group=TRUE,console=TRUE)
write.table(out_CE$groups, file = "out_CE.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_CO<-aov(CO~groups, data=fisicoq); summary(model_CO)
out_CO <- HSD.test(model_CO,"groups", group=TRUE,console=TRUE)
write.table(out_CO$groups, file = "out_CO.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_CT<-aov(CT~groups, data=fisicoq); summary(model_CT)
out_CT <- HSD.test(model_CT,"groups", group=TRUE,console=TRUE)
write.table(out_CT$groups, file = "out_CT.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_NT<-aov(NT~groups, data=fisicoq); summary(model_NT)
out_NT <- HSD.test(model_NT,"groups", group=TRUE,console=TRUE)
write.table(out_NT$groups, file = "out_NT.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Ca<-aov(Ca~groups, data=fisicoq); summary(model_Ca)
out_Ca <- HSD.test(model_Ca,"groups", group=TRUE,console=TRUE)
write.table(out_Ca$groups, file = "out_Ca.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_K<-aov(K~groups, data=fisicoq); summary(model_K)
out_K <- HSD.test(model_K,"groups", group=TRUE,console=TRUE)
write.table(out_K$groups, file = "out_K.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Mg<-aov(Mg~groups, data=fisicoq); summary(model_Mg)
out_Mg <- HSD.test(model_Mg,"groups", group=TRUE,console=TRUE)
write.table(out_Mg$groups, file = "out_Mg.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Na<-aov(Na~groups, data=fisicoq); summary(model_Na)
out_Na <- HSD.test(model_Na,"groups", group=TRUE,console=TRUE)
write.table(out_Na$groups, file = "out_Na.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_CICE<-aov(CICE~groups, data=fisicoq); summary(model_CICE)
out_CICE <- HSD.test(model_CICE,"groups", group=TRUE,console=TRUE)
write.table(out_CICE$groups, file = "out_CICE.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_P<-aov(P~groups, data=fisicoq); summary(model_P)
out_P <- HSD.test(model_P,"groups", group=TRUE,console=TRUE)
write.table(out_P$groups, file = "out_P.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_S<-aov(S~groups, data=fisicoq); summary(model_S)
out_S <- HSD.test(model_S,"groups", group=TRUE,console=TRUE)
write.table(out_S$groups, file = "out_S.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Cu<-aov(Cu~groups, data=fisicoq); summary(model_Cu)
out_Cu <- HSD.test(model_Cu,"groups", group=TRUE,console=TRUE)
write.table(out_Cu$groups, file = "out_Cu.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Fe<-aov(Fe~groups, data=fisicoq); summary(model_Fe)
out_Fe <- HSD.test(model_Fe,"groups", group=TRUE,console=TRUE)
write.table(out_Fe$groups, file = "out_Fe.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Mn<-aov(Mn~groups, data=fisicoq); summary(model_Mn)
out_Mn <- HSD.test(model_Mn,"groups", group=TRUE,console=TRUE)
write.table(out_Mn$groups, file = "out_Mn.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Zn<-aov(Zn~groups, data=fisicoq); summary(model_Zn)
out_Zn <- HSD.test(model_Zn,"groups", group=TRUE,console=TRUE)
write.table(out_Zn$groups, file = "out_Zn.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_B<-aov(B~groups, data=fisicoq); summary(model_B)
out_B <- HSD.test(model_B,"groups", group=TRUE, console=TRUE)
write.table(out_B$groups, file = "out_B.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Ar<-aov(Arcilla~groups, data=fisicoq); summary(model_Ar)
out_Ar <- HSD.test(model_Ar,"groups", group=TRUE, console=TRUE)
write.table(out_Ar$groups, file = "out_Ar.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_L <-aov(Limo~groups, data=fisicoq); summary(model_L)
out_L <- HSD.test(model_L,"groups", group=TRUE, console=TRUE)
write.table(out_L$groups, file = "out_L.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_A <-aov(Arena~groups, data=fisicoq); summary(model_A)
out_A <- HSD.test(model_A,"groups", group=TRUE, console=TRUE)
write.table(out_A$groups, file = "out_A.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_NO3 <-aov(NO3~groups, data=fisicoq); summary(model_A)
out_NO3 <- HSD.test(model_NO3,"groups", group=TRUE, console=TRUE)
write.table(out_NO3$groups, file = "out_NO3.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_NH4 <-aov(NH4~groups, data=fisicoq); summary(model_A)
out_NH4 <- HSD.test(model_NH4,"groups", group=TRUE, console=TRUE)
write.table(out_NH4$groups, file = "out_NH4.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

getwd()  # print the current working directory
setwd("~RamosFigueroa/Documents/OneDrive/Work/PROYECTO MANGLAR/BIOINF_ANALYSIS_MANGROVE_16S/Segundo Muestreo 16S/CanjerejosCamarones")
