library(plotrix)
library(agricolae)

fisicoq <- read.delim("/Users/RamosFigueroa/Documents/OneDrive/Work/Tesis\ Maestr??a\ Microbiolog??a/Documentos\ de\ Tesis/fisicoquimicos\ 2015.txt")
fisicoq_promedios <- round(aggregate(fisicoq[, 3:23], list(fisicoq$group), mean), 2)
fisicoq_errorstand <- round(aggregate(fisicoq[, 3:23], list(fisicoq$group), std.error), 2)
fisicoq_promediosT <- t(fisicoq_promedios)
fisicoq_errorstandT <- t(fisicoq_errorstand)

write.table(fisicoq_promedios, file = "/Users/RamosFigueroa/Documents/OneDrive/Work/Tesis\ Maestr??a\ Microbiolog??a/Documentos\ de\ Tesis/fisicoq_promedios.txt", sep = "\t", quote = FALSE, col.names = TRUE)
write.table(fisicoq_errorstand, file = "/Users/RamosFigueroa/Documents/OneDrive/Work/Tesis\ Maestr??a\ Microbiolog??a/Documentos\ de\ Tesis/fisicoq_errorstand.txt", sep = "\t", quote = FALSE, col.names = TRUE)

write.table(fisicoq_promediosT, file = "/Users/RamosFigueroa/Documents/OneDrive/Work/Tesis\ Maestr??a\ Microbiolog??a/Documentos\ de\ Tesis/fisicoq_promediosT.txt", sep = "\t", quote = FALSE, col.names = TRUE)
write.table(fisicoq_errorstandT, file = "/Users/RamosFigueroa/Documents/OneDrive/Work/Tesis\ Maestr??a\ Microbiolog??a/Documentos\ de\ Tesis/fisicoq_errorstandT.txt", sep = "\t", quote = FALSE, col.names = TRUE)


# ANOVAS Y TUKEY

model_pH<-aov(pH~group, data=fisicoq); summary(model_pH)
out_pH <- HSD.test(model_pH,"group", group=TRUE,console=TRUE)

model_CE<-aov(CE~group, data=fisicoq); summary(model_CE)
out_CE <- HSD.test(model_CE,"group", group=TRUE,console=TRUE)

model_CO<-aov(CO~group, data=fisicoq); summary(model_CO)
out_CO <- HSD.test(model_CO,"group", group=TRUE,console=TRUE)

model_CT<-aov(CT~group, data=fisicoq); summary(model_CT)
out_CT <- HSD.test(model_CT,"group", group=TRUE,console=TRUE)

model_NT<-aov(NT~group, data=fisicoq); summary(model_NT)
out_NT <- HSD.test(model_NT,"group", group=TRUE,console=TRUE)

model_Ca<-aov(Ca~group, data=fisicoq); summary(model_Ca)
out_Ca <- HSD.test(model_Ca,"group", group=TRUE,console=TRUE)

model_K<-aov(K~group, data=fisicoq); summary(model_K)
out_K <- HSD.test(model_K,"group", group=TRUE,console=TRUE)

model_Mg<-aov(Mg~group, data=fisicoq); summary(model_Mg)
out_Mg <- HSD.test(model_Mg,"group", group=TRUE,console=TRUE)

model_Na<-aov(Na~group, data=fisicoq); summary(model_Na)
out_Na <- HSD.test(model_Na,"group", group=TRUE,console=TRUE)

model_CICE<-aov(CICE~group, data=fisicoq); summary(model_CICE)
out_CICE <- HSD.test(model_CICE,"group", group=TRUE,console=TRUE)

model_P<-aov(P~group, data=fisicoq); summary(model_P)
out_P <- HSD.test(model_P,"group", group=TRUE,console=TRUE)

model_S<-aov(S~group, data=fisicoq); summary(model_S)
out_S <- HSD.test(model_S,"group", group=TRUE,console=TRUE)

model_Cu<-aov(Cu~group, data=fisicoq); summary(model_Cu)
out_Cu <- HSD.test(model_Cu,"group", group=TRUE,console=TRUE)

model_Fe<-aov(Fe~group, data=fisicoq); summary(model_Fe)
out_Fe <- HSD.test(model_Fe,"group", group=TRUE,console=TRUE)

model_Mn<-aov(Mn~group, data=fisicoq); summary(model_Mn)
out_Mn <- HSD.test(model_Mn,"group", group=TRUE,console=TRUE)

model_Zn<-aov(Zn~group, data=fisicoq); summary(model_Zn)
out_Zn <- HSD.test(model_Zn,"group", group=TRUE,console=TRUE)

model_B<-aov(B~group, data=fisicoq); summary(model_B)
out_B <- HSD.test(model_B,"group", group=TRUE, console=TRUE)

model_Ar<-aov(Ar~group, data=fisicoq); summary(model_Ar)
out_Ar <- HSD.test(model_Ar,"group", group=TRUE, console=TRUE)

model_L<-aov(L~group, data=fisicoq); summary(model_L)
out_L <- HSD.test(model_L,"group", group=TRUE, console=TRUE)

model_A<-aov(A~group, data=fisicoq); summary(model_A)
out_A <- HSD.test(model_A,"group", group=TRUE, console=TRUE)

