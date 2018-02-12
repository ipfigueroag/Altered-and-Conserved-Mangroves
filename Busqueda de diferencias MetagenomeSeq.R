######## MetagenomeSeq Segundo Muestreo 16S #######
####### Ingrid Paola's Customized Scripts ########

library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(plotly)
library(metagenomeSeq)

getwd()  # print the current working directory
setwd("~RamosFigueroa/Documents/OneDrive/Work/PROYECTO MANGLAR/BIOINF_ANALYSIS_MANGROVE_16S/Segundo Muestreo 16S/CanjerejosCamarones")

otu <- "cp_l_0.1_otu_counts.norm.rds";   #subiendo tabla OTUs
otu_tab <- readRDS(otu); otu_tab <- as.data.frame(otu_tab)   #leyendo archivo .rds y conversion en dataframe
write.table(otu_tab, file="tabla_otus.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE) #exportando datos desde R a ubicaci??n local

tax <- "cp_l_0.1_taxonomy.rds"  #subiendo tabla taxonomy
tax_tab <- readRDS(tax); #convertir objeto en tabla
write.table(tax_tab, file="tabla_taxonomia.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE) #exportando datos desde R a ubicaci??n local

tax_tab <- tax_tab[-c(1:2,9:14)]  # eliminar columnas

metadata <- "metadata.txt"; metadata <- read.delim(metadata) #subiendo metadata
metadata <- sample_data(metadata); 

### Load Data to MetagenoSeq ####

dataDirectory <- list.files(path = "~RamosFigueroa/Documents/OneDrive/Work/PROYECTO MANGLAR/BIOINF_ANALYSIS_MANGROVE_16S/Segundo Muestreo 16S/CanjerejosCamarones")
otu <- load_meta("tabla_otus.txt")
taxa <- read.delim("tabla_taxonomia.txt")
ambiental <- load_phenoData("metadata.txt")

phenotypeData <- AnnotatedDataFrame(ambiental); phenotypeData
OTUdata <- AnnotatedDataFrame(taxa); OTUdata
obj <- newMRexperiment(otu$counts,phenoData=phenotypeData,featureData=OTUdata) #Se requiere para an??lisis posteriores


##### Permutaciones y graficos CAP ###########

#a_group <- adonis(dist.bc ~ groups, data = metadata, perm = 9999); a_CE <- adonis(dist.bc ~ Cond, data = metadata, perm = 9999); a_CO <- adonis(dist.bc ~ Corg, data = metadata, perm = 9999); a_pH <- adonis(dist.bc ~ pH, data = metadata, perm = 9999); a_NT <- adonis(dist.bc ~ NT, data = metadata, perm = 9999);
#a_group; a_CE; a_CO; a_Arcilla; a_pH; a_NT

metadata <- sample_data(metadata)
row.names(metadata) <- metadata$ID
obphyloeq_merge <- merge_phyloseq(obphyloseq, metadata)

GP.ord <- ordinate(obphyloeq_merge, "DCA", "bray")
p2 = plot_ordination(obphyloeq_merge, GP.ord, type = "samples", colo = "groups", title="Sample"); 
p1 = plot_ordination(obphyloeq_merge, GP.ord, color="groups", title="taxa")
print(p1); print(p2)

## Constrained Ordinations ##

# Remove data points with missing metadata
moth_not_na <- obphyloeq_merge %>%
  subset_samples(
      !is.na(Plant_rizosfera) &
      !is.na(CE) 
  )

bray_not_na <- phyloseq::distance(physeq = moth_not_na, method = "bray")


# CAP ordinate
cap_ord <- ordinate(
  physeq = moth_not_na, 
  method = "CAP",
  distance = bray_not_na,
  formula = ~ Plant_rizosfera + CE
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = moth_not_na, 
  ordination = cap_ord,
  color = NULL,
  label =  "ID",
  axes = c(1,2)
) + 
  aes(shape = Station) + 
  geom_point(aes(colour = NULL), alpha = 0.4, size = 15) + 
  geom_point(colour = "black", size = 2) + theme(labels(element_text(size = 16))) + 
scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a", 
"#1919ff", "darkorchid3", "magenta", "blue", "red", "#5E738F","#D1A33D", "#8A7C64", "#599861", "4", "3", "#CBD588", "#5F7FC7", "orange", "#DA5724", "#508578", "#CD9BCD", "#AD6F3B", "#673770","#D14285", "#C84248", "#652926", "#8569D5", "8", "#4569D2", "11", "22")
)


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = 1, 
    data = arrowdf, 
    color = "black", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 6,  
    data = arrowdf, 
    show.legend = FALSE
  ) + theme(axis.text=element_text(size=15),           #tama??o del texto de ejes
            axis.title=element_text(size=18)) + theme(legend.text=element_text(size=16))

anova(cap_ord)

# Black & White #
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = 1, 
    data = arrowdf, 
    color = "black", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 15,  
    data = arrowdf, 
    show.legend = TRUE
  ) + theme(axis.text=element_text(size=25),           #tama??o del texto de ejes
            axis.title=element_text(size=30))


### Diferencias Significativas en Phylum y Generos ##

library(plotrix)
library(agricolae)

getwd()  # print the current working directory
setwd("~RamosFigueroa/Documents/OneDrive/Work/PROYECTO\ MANGLAR/BIOINF_ANALYSIS_MANGROVE_16S/Segundo\ Muestreo\ 16S/CanjerejosCamarones/Resultados\ Paola\ Figueroa/")



generos <- read.delim("Transpuesta Generos Organizados con diferencias sig STAMP.txt"); generos <- as.data.frame(generos)

generos_promedios1 <- aggregate(generos[, 3:57], list(generos$Groups), mean); row.names(generos_promedios1) <- generos_promedios1[,1]
generos_promedios <- round(generos_promedios1[, 2:56], 2)

generos_errorstand1 <- aggregate(generos[, 3:57], list(generos$Groups), std.error); row.names(generos_errorstand1) <- generos_errorstand1[,1]
generos_errorstand <- round(generos_errorstand1[, 2:56], 2)



write.table(generos_promedios, file = "generos_promedios.txt", sep = "\t", quote = FALSE, col.names = TRUE)
write.table(generos_errorstand, file = "generos_errorstand.txt", sep = "\t", quote = FALSE, col.names = TRUE)

#write.table(fisicoq_promediosT, file = "/Users/RamosFigueroa/Documents/OneDrive/Work/Tesis\ Maestr??a\ Microbiolog??a/Documentos\ de\ Tesis/fisicoq_promediosT.txt", sep = "\t", quote = FALSE, col.names = TRUE)
#write.table(fisicoq_errorstandT, file = "/Users/RamosFigueroa/Documents/OneDrive/Work/Tesis\ Maestr??a\ Microbiolog??a/Documentos\ de\ Tesis/fisicoq_errorstandT.txt", sep = "\t", quote = FALSE, col.names = TRUE)


# ANOVAS Y TUKEY

model_Halobacillus <-aov(Halobacillus~Groups, data=generos); summary(model_Halobacillus)
out_Halobacillus <- HSD.test(model_Halobacillus,"Groups", group=TRUE,console=TRUE)
write.table(out_Halobacillus$groups, file = "out_Halobacillus.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Salinibacter <-aov(Salinibacter~Groups, data=generos); summary(model_Salinibacter)
out_Salinibacter <- HSD.test(model_Salinibacter,"Groups", group=TRUE,console=TRUE)
write.table(out_Salinibacter$groups, file = "out_Salinibacter.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Pir4_lineage <-aov(Pir4_lineage~Groups, data=generos); summary(model_Pir4_lineage)
out_Pir4_lineage <- HSD.test(model_Pir4_lineage,"Groups", group=TRUE,console=TRUE)
write.table(out_Pir4_lineage$groups, file = "out_Pir4_lineage.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Pontibacillus <-aov(Pontibacillus~Groups, data=generos); summary(model_Pontibacillus)
out_Pontibacillus <- HSD.test(model_Pontibacillus,"Groups", group=TRUE,console=TRUE)
write.table(out_Pontibacillus$groups, file = "out_Pontibacillus.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Bythopirellula <-aov(Bythopirellula~Groups, data=generos); summary(model_Bythopirellula)
out_Bythopirellula <- HSD.test(model_Bythopirellula,"Groups", group=TRUE,console=TRUE)
write.table(out_Bythopirellula$groups, file = "out_Bythopirellula.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Clostridium_sensu_stricto <-aov(Clostridium_sensu_stricto~Groups, data=generos); summary(model_Clostridium_sensu_stricto)
out_Clostridium_sensu_stricto <- HSD.test(model_Clostridium_sensu_stricto,"Groups", group=TRUE,console=TRUE)
write.table(out_Clostridium_sensu_stricto$groups, file = "out_Clostridium_sensu_stricto.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Nitrolancea <-aov(Nitrolancea~Groups, data=generos); summary(model_Nitrolancea)
out_Nitrolancea <- HSD.test(model_Nitrolancea,"Groups", group=TRUE,console=TRUE)
write.table(out_Nitrolancea$groups, file = "out_Nitrolancea.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Dichotomicrobium<-aov(Dichotomicrobium~Groups, data=generos); summary(model_Dichotomicrobium)
out_Dichotomicrobium <- HSD.test(model_Dichotomicrobium,"Groups", group=TRUE,console=TRUE)
write.table(out_Dichotomicrobium$groups, file = "out_Dichotomicrobium.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Sva0081_sediment_group<-aov(Sva0081_sediment_group~Groups, data=generos); summary(model_Sva0081_sediment_group)
out_Sva0081_sediment_group <- HSD.test(model_Sva0081_sediment_group,"Groups", group=TRUE,console=TRUE)
write.table(out_Sva0081_sediment_group$groups, file = "out_Sva0081_sediment_group.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Urania.1B.19_marine_sediment_group<-aov(Urania.1B.19_marine_sediment_group~Groups, data=generos); summary(model_Urania.1B.19_marine_sediment_group)
out_Urania.1B.19_marine_sediment_group <- HSD.test(model_Urania.1B.19_marine_sediment_group,"Groups", group=TRUE,console=TRUE)
write.table(out_Urania.1B.19_marine_sediment_group$groups, file = "out_Urania.1B.19_marine_sediment_group.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Blastopirellula<-aov(Blastopirellula~Groups, data=generos); summary(model_Blastopirellula)
out_Blastopirellula <- HSD.test(model_Blastopirellula,"Groups", group=TRUE,console=TRUE)
write.table(out_Blastopirellula$groups, file = "out_Blastopirellula.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Candidatus_Thiobios<-aov(Candidatus_Thiobios~Groups, data=generos); summary(model_Candidatus_Thiobios)
out_Candidatus_Thiobios <- HSD.test(model_Candidatus_Thiobios,"Groups", group=TRUE,console=TRUE)
write.table(out_Candidatus_Thiobios$groups, file = "out_Candidatus_Thiobios.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Methyloceanibacter<-aov(Methyloceanibacter~Groups, data=generos); summary(model_Methyloceanibacter)
out_Methyloceanibacter <- HSD.test(model_Methyloceanibacter,"Groups", group=TRUE,console=TRUE)
write.table(out_Methyloceanibacter$groups, file = "out_Methyloceanibacter.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Paenibacillus<-aov(Paenibacillus~Groups, data=generos); summary(model_Paenibacillus)
out_Paenibacillus <- HSD.test(model_Paenibacillus,"Groups", group=TRUE,console=TRUE)
write.table(out_Paenibacillus$groups, file = "out_Paenibacillus.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Crossiella <-aov(Crossiella~Groups, data=generos); summary(model_Crossiella)
out_Crossiella <- HSD.test(model_Crossiella,"Groups", group=TRUE,console=TRUE)
write.table(out_Crossiella$groups, file = "out_Crossiella.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Nitrosococcus <-aov(Nitrosococcus~Groups, data=generos); summary(model_Nitrosococcus)
out_Nitrosococcus <- HSD.test(model_Nitrosococcus,"Groups", group=TRUE,console=TRUE)
write.table(out_Nitrosococcus$groups, file = "out_Nitrosococcus.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Portibacter<-aov(Portibacter~Groups, data=generos); summary(model_Portibacter)
out_Portibacter <- HSD.test(model_Portibacter,"Groups", group=TRUE, console=TRUE)
write.table(out_Portibacter$groups, file = "out_Portibacter.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Tumebacillus<-aov(Tumebacillus~Groups, data=generos); summary(model_Tumebacillus)
out_Tumebacillus <- HSD.test(model_Tumebacillus,"Groups", group=TRUE, console=TRUE)
write.table(out_Tumebacillus$groups, file = "out_Tumebacillus.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Rhodopirellula <-aov(Rhodopirellula~Groups, data=generos); summary(model_Rhodopirellula)
out_Rhodopirellula <- HSD.test(model_Rhodopirellula,"Groups", group=TRUE, console=TRUE)
write.table(out_Rhodopirellula$groups, file = "out_Rhodopirellula.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Deferrisoma<-aov(Deferrisoma~Groups, data=generos); summary(model_Deferrisoma)
out_Deferrisoma <- HSD.test(model_Deferrisoma,"Groups", group=TRUE, console=TRUE)
write.table(out_Deferrisoma$groups, file = "out_Deferrisoma.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Oceanobacillus<-aov(Oceanobacillus~Groups, data=generos); summary(model_Oceanobacillus)
out_Oceanobacillus <- HSD.test(model_Oceanobacillus,"Groups", group=TRUE, console=TRUE)
write.table(out_Oceanobacillus$groups, file = "out_Oceanobacillus.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

model_Rhodomicrobium<-aov(Rhodomicrobium~Groups, data=generos); summary(model_Rhodomicrobium)
out_Rhodomicrobium <- HSD.test(model_Rhodomicrobium,"Groups", group=TRUE, console=TRUE)
write.table(out_Rhodomicrobium$groups, file = "out_Rhodomicrobium.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

