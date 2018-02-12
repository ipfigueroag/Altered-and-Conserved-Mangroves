
### phyloseq ####
### Customized Scripts of Ingrid Paola Figueroa Galvis #####

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
opt <- "cp_l_0.1_otu_counts.norm.rds"  #subiendo tabla OTUs
df.r <- readRDS(opt); df <- as.data.frame(df.r) #convertir objeto en tabla y luego transformar a dataframe
prom1 <- ((df$`1`+df$`4`+df$`8`)/3); prom2 <- ((df$`10`+df$`13`+df$`16`)/3); prom3 <- ((df$`17`+df$`18`+df$`19`)/3);
          prom4 <- ((df$`20`+df$`21`+df$`22`)/3); prom5 <- ((df$`23`+df$`24`+df$`25`)/3); prom6 <- ((df$`26`+df$`27`+df$`28`)/3);
tab_prom <- data.frame(`VC1`=prom1, `VC2`=prom2, `VC3`=prom3, `PNC1`=prom4, `PNC2`=prom5, `PNC3`=prom6) #creando tabla
row.names(tab_prom)<- row.names(df)  #colocando otus en rownames de nueva tabla

obtax <- "cp_l_0.1_taxonomy.rds"  #subiendo tabla taxonomy
ttax <- readRDS(obtax); ttax <- ttax[-c(1:2,9:14)]; ttax <- as.matrix(ttax)   #convertir objeto en tabla y luego transformar a dataframe

ph <- "cp_l_0.1_phenotype.rds"; phen <- readRDS(ph) #subiendo metadata

OTU = otu_table(df, taxa_are_rows = TRUE)
TAX = tax_table(ttax); row.names(TAX) <- row.names(ttax); colnames(TAX) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

obphyloseq <- phyloseq(OTU, TAX) # construyendo objeto clase phyloseq desde componentes
phen <- sample_data(phen) # cargando medatada a phyloseq

phylum_colors <- c("#5E738F","#D1A33D", "#8A7C64", "#599861", "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD", "#AD6F3B", "#673770","#D14285", "#C84248", "#652926", "#8569D5", "8", "#4569D2", "11", "22", "4", "3", "9", "10", "6", "#5E738F","#D1A33D", "#8A7C64", "#599861", "#CBD588", "#5F7FC7", "orange","#DA5724", "#5E738F","#D1A33D", "#8A7C64", "#599861", "#CBD588", "#5F7FC7", "orange","#DA5724")

######### Ploting Phylum 16S #########

phylum2 <- obphyloseq  %>% 
  tax_glom(taxrank = "Phylum")  %>%                  # agglomerate at phylum level
  #transform_sample_counts(function(x) {100 *x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.5) %>%                         # Filter out low abundance taxa
  arrange(Phylum) 
View(phylum2)

write.table(phylum2, file="Conteos_Phylum.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE) #exportando datos desde R a ubicaci??n local

plotbars <- ggplot(phylum2, aes(x = Sample, y = Abundance, fill = Phylum, order = sort(phylum2$Phylum, decreasing = TRUE))) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = phylum_colors) +
  theme(axis.text=element_text(size=40, colour = "black", family = "Times"),           #tama??o del texto de ejes
        axis.title=element_text(size=50, family = "Times"),
        legend.text=element_text(size=50, face = "italic", family = "Times")) +        #tama??o t??tulo eje x
  scale_x_discrete(labels=c("NP1", "NP2", "NP3", "VC1", "VC2", "VC3")) +
  scale_y_continuous(labels = percent_format()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 3, keyheight = 3)) +
  ylab("\nRelative Abundance > 0.5%") + 
  xlab("\nSamples\n")
 
plotbars

### Class 16S #####

Class16S <- obphyloseq  %>% 
  tax_glom(taxrank = "Class")  %>%                  # agglomerate at phylum level
  #transform_sample_counts(function(x) {100 *x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 1.5) %>%                         # Filter out low abundance taxa
  arrange(Class) 
View(Class16S)

write.table(Class16S, file="Conteos_Class16S.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE) #exportando datos desde R a ubicaci??n local

plotbarsClass <- ggplot(Class16S, aes(x = Sample, y = Abundance, fill = Class, order = sort(Class16S$Class, decreasing = TRUE))) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = phylum_colors) +
  theme(axis.text=element_text(size=30, colour = "black", family = "Times"),           #tama??o del texto de ejes
        axis.title=element_text(size=30, family = "Times"),
        legend.text=element_text(size=30, face = "italic", family = "Times")) +        #tama??o t??tulo eje x
  scale_x_discrete(labels=c("NP1", "NP2", "NP3", "VC1", "VC2", "VC3")) +
  scale_y_continuous(labels = percent_format()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 2, keyheight = 2)) +
  ylab("\nRelative Abundance > 1.5%") + 
  xlab("\nSamples\n")

plotbarsClass


### Order 16S #####

Order16S <- obphyloseq  %>% 
  tax_glom(taxrank = "Order")  %>%                  # agglomerate at phylum level
  #transform_sample_counts(function(x) {100 *x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 2.5) %>%                         # Filter out low abundance taxa
  arrange(Order) 
View(Order16S)

write.table(Order16S, file="Conteos_Order16S.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE) #exportando datos desde R a ubicaci??n local

plotbarsOrder <- ggplot(Order16S, aes(x = Sample, y = Abundance, fill = Order, order = sort(Order16S$Order, decreasing = TRUE))) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = phylum_colors) +
  theme(axis.text=element_text(size=30, colour = "black", family = "Times"),           #tama??o del texto de ejes
        axis.title=element_text(size=30, family = "Times"),
        legend.text=element_text(size=30, face = "italic", family = "Times")) +        #tama??o t??tulo eje x
  scale_x_discrete(labels=c("14.61", "2.80", "61.52")) +
  scale_y_continuous(labels = percent_format()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 2, keyheight = 2)) +
  ylab("\nRelative Abundance > 2.5%") + 
  xlab("\nSalinity (ppt)")

plotbarsOrder


### Family 16S #####

Family16S <- obphyloseq  %>% 
  tax_glom(taxrank = "Family")  %>%                  # agglomerate at phylum level
  #transform_sample_counts(function(x) {100 *x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 2.5) %>%                         # Filter out low abundance taxa
  arrange(Family) 
View(Family16S)

write.table(Family16S, file="Conteos_family16S.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE) #exportando datos desde R a ubicaci??n local

plotbarsFamily <- ggplot(Family16S, aes(x = Sample, y = Abundance, fill = Family, order = sort(Family16S$Family, decreasing = TRUE))) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = phylum_colors) +
  theme(axis.text=element_text(size=30, colour = "black", family = "Times"),           #tama??o del texto de ejes
        axis.title=element_text(size=30, family = "Times"),
        legend.text=element_text(size=30, face = "italic", family = "Times")) +        #tama??o t??tulo eje x
  scale_x_discrete(labels=c("14.61", "2.80", "61.52")) +
  scale_y_continuous(labels = percent_format()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 2, keyheight = 2)) +
  ylab("\nRelative Abundance > 2.5%") + 
  xlab("\nSalinity (ppt)")

plotbarsFamily


### Genus 16S #####

Genus16S <- obphyloseq  %>% 
  tax_glom(taxrank = "Genus")  %>%                  # agglomerate at phylum level
  #transform_sample_counts(function(x) {100 *x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 2) %>%                         # Filter out low abundance taxa
  arrange(Genus) 
View(Genus16S)

write.table(Genus16S, file="Conteos_Genus16S.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE) #exportando datos desde R a ubicaci??n local

plotbarsGenus <- ggplot(Genus16S, aes(x = Sample, y = Abundance, fill = Genus, order = sort(Genus16S$Genus, decreasing = TRUE))) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = phylum_colors) +
  theme(axis.text=element_text(size=30, colour = "black", family = "Times"),           #tama??o del texto de ejes
        axis.title=element_text(size=30, family = "Times"),
        legend.text=element_text(size=25, face = "italic", family = "Times")) +        #tama??o t??tulo eje x
  scale_x_discrete(labels=c("14.61", "2.80", "61.52")) +
  scale_y_continuous(labels = percent_format()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 2, keyheight = 2)) +
  ylab("\nRelative Abundance > 2%") + 
  xlab("\nSalinity (ppt)")

plotbarsGenus

###### rarefaction curve #####

library(RColorBrewer)

df_round <- round(df, digits = 0); #redondeando los valores del dataframe de OTUs a numeros enteros
metadato <- read.table("metadata.txt", header = TRUE) #subiendo metadata

#OTU2 = otu_table(df_round, taxa_are_rows = TRUE)
#objphyloseq2 <- phyloseq(OTU2, TAX) # construyendo objeto clase phyloseq desde componentes
#met <- sample_data(metadato) # cargando medatada a phyloseq
#amp_rarecurve(objphyloseq2)  # creando curva de rarefaccion
 
metadato <- sample_data(metadato);
rownames(metadato) <- metadato$ID; # Assign rownames to be Sample ID's
objphyloseq2 <- merge_phyloseq(objphyloseq2, metadato) # Merge mothurdata object with sample metadata

sample_data(objphyloseq2)$groups <- c("VC1", "VC2", "VC2", "VC2", "VC3", "VC3", "VC3", "NP1", "NP1", "NP1", "NP2", "NP2", "NP2", "NP3", "NP3", "NP3", "VC1", "VC1")
sample_data(objphyloseq2)$groups

p <- ggrare(objphyloseq2, step = 100, color = "groups", se = FALSE)
plot(p + theme_bw(base_size = 40, base_family = "Times") +
     theme(axis.text=element_text(colour = "black"), legend.key.size = unit(0.8,"cm"), 
           legend.key.height = unit(1.2,"cm"), legend.title = element_text(vjust = 0,5) ) +
     geom_line(size=2) +
  ylab("OTUs") + 
  xlab("Reads Sampled"))
    


#plotbars <- ggplot(p, aes(x = Sample, y = Abundance, fill = Phylum, order = sort(phylum2$Phylum, decreasing = TRUE))) + 
 # geom_bar(stat = "identity", position = "fill")

#### Guardar datos de conteos en diferentes niveles ######

write.table(ttax, file="countstax.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE) #exportando datos desde R a ubicaci??n local

write.table(df, file="countsotus.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE) #exportando datos desde R a ubicaci??n local

write.table(phen, file="metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE) #exportando datos desde R a ubicaci??n local


##### Congreso Biotecnologia ##### METAGENOMA ... phylum #####
setwd("~RamosFigueroa/Documents/OneDrive/Work/PROYECTO\ MANGLAR/BIOINF_ANALYSIS_MANGROVE_16S/Comparison\ Metagenome\ AND\ PICRUSt/")

taxmeta <- "phylum_tax.txt"; taxmeta <- read.delim(taxmeta)  #subiendo tabla OTUs
taxmeta <- as.data.frame(taxmeta) #convertir objeto en tabla y luego transformar a dataframe
prom1meta <- ((taxmeta$gi_4A_R1_80+taxmeta$gi_4A_R2_80+taxmeta$gi_4B_R1_80+taxmeta$gi_4B_R2_80+taxmeta$gi_4C_R1_80+taxmeta$gi_4C_R2_80)/6); prom2meta <- ((taxmeta$gi_2A_R1_80+taxmeta$gi_2A_R2_80+taxmeta$gi_2B_R1_80+taxmeta$gi_2B_R2_80+taxmeta$gi_2C_R1_80+taxmeta$gi_2C_R2_80)/6); prom3meta <- ((taxmeta$gi_3A_R1_80+taxmeta$gi_3A_R2_80+taxmeta$gi_3B_R1_80+taxmeta$gi_3B_R2_80+taxmeta$gi_3C_R1_80+taxmeta$gi_3C_R2_80)/6)
tab_prommeta <- data.frame(Site2=prom2meta, Site3=prom3meta, Site4=prom1meta) #creando tabla
#tabla taxonomy:
otumeta <- data.frame(Phylum=taxmeta$X.Datasets); row.names(otumeta) <- c("Otu1", "Otu2", "Otu3", "Otu4", "Otu5", "Otu6", "Otu7", "Otu8", "Otu9", "Otu10", "Otu11", "Otu12", "Otu13", "Otu14", "Otu15", "Otu16", "Otu17", "Otu18", "Otu19", "Otu20", "Otu21", "Otu22", "Otu23", "Otu24", "Otu25", "Otu26", "Otu27", "Otu28", "Otu29", "Otu30", "Otu31", "Otu32", "Otu33", "Otu34", "Otu35", "Otu36", "Otu37")   #creando tabla de OTUs y colocando valores en rownames
row.names(tab_prommeta)<- row.names(otumeta)  #colocando otus en rownames tabla de otus

## Cargando objeto de phyloseq
OTUMETAG = otu_table(tab_prommeta, taxa_are_rows = TRUE)
TAXMETAG = tax_table(otumeta); row.names(TAXMETAG) <- row.names(otumeta); colnames(TAXMETAG) <- c("Phylum")
obphyloseqMETAG <- phyloseq(OTUMETAG, TAXMETAG)

phylumMETAG <- obphyloseqMETAG  %>% 
  tax_glom(taxrank = "Phylum")  %>%                  # agglomerate at phylum level
  transform_sample_counts(function(x) {100 *x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.5) %>%                         # Filter out low abundance taxa
  arrange(Phylum)
View(phylumMETAG)

write.table(phylumMETAG, file="RelativeAbundancePhylumMETAGENOMA.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE) #exportando datos desde R a ubicaci??n local

plotbarsMETAG <- ggplot(phylumMETAG, aes(x = Sample, y = Abundance, fill = Phylum, order = sort(phylumMETAG$Phylum, decreasing = TRUE))) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = phylum_colors) +
  theme(axis.text=element_text(size=15, colour = "black", family = "Times"),           #tama??o del texto de ejes
        axis.title=element_text(size=15, family = "Times"),
        legend.text=element_text(size=15, face = "italic", family = "Times")) +        #tama??o t??tulo eje x
  scale_y_continuous(labels = percent_format()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("\nRelative Abundance") + 
  xlab("\nSites")

plotbarsMETAG


##### Congreso Biotecnologia ##### METAGENOMA ... class #####

taxmetaclass <- "class_tax.txt"; taxmetaclass <- read.delim(taxmetaclass)  #subiendo tabla OTUs
taxmetaclass <- as.data.frame(taxmetaclass) #convertir objeto en tabla y luego transformar a dataframe
prom1metaclass <- ((taxmetaclass$gi_4A_R1_80+taxmetaclass$gi_4A_R2_80+taxmetaclass$gi_4B_R1_80+taxmetaclass$gi_4B_R2_80+taxmetaclass$gi_4C_R1_80+taxmetaclass$gi_4C_R2_80)/6); prom2metaclass <- ((taxmetaclass$gi_2A_R1_80+taxmetaclass$gi_2A_R2_80+taxmetaclass$gi_2B_R1_80+taxmetaclass$gi_2B_R2_80+taxmetaclass$gi_2C_R1_80+taxmetaclass$gi_2C_R2_80)/6); prom3metaclass <- ((taxmetaclass$gi_3A_R1_80+taxmetaclass$gi_3A_R2_80+taxmetaclass$gi_3B_R1_80+taxmetaclass$gi_3B_R2_80+taxmetaclass$gi_3C_R1_80+taxmetaclass$gi_3C_R2_80)/6)
tab_prommetaclass <- data.frame(Site2=prom2metaclass, Site3=prom3metaclass, Site4=prom1metaclass) #creando tabla
#tabla taxonomy:
otumetaclass <- data.frame(Class=taxmetaclass$X.Datasets); row.names(otumetaclass) <- c("Otu1", "Otu2", "Otu3", "Otu4", "Otu5", "Otu6", "Otu7", "Otu8", "Otu9", "Otu10", "Otu11", "Otu12", "Otu13", "Otu14", "Otu15", "Otu16", "Otu17", "Otu18", "Otu19", "Otu20", "Otu21", "Otu22", "Otu23", "Otu24", "Otu25", "Otu26", "Otu27", "Otu28", "Otu29", "Otu30", "Otu31", "Otu32", "Otu33", "Otu34", "Otu35", "Otu36", "Otu37", "Otu38", "Otu39", "Otu40", "Otu41", "Otu42", "Otu43", "Otu44", "Otu45", "Otu46", "Otu47", "Otu48", "Otu49", "Otu50", "Otu51", "Otu52", "Otu53", "Otu54", "Otu55", "Otu56", "Otu57", "Otu58", "Otu59", "Otu60", "Otu61", "Otu62", "Otu63", "Otu64", "Otu65", "Otu66", "Otu67", "Otu68", "Otu69", "Otu70", "Otu71")   #creando tabla de OTUs y colocando valores en rownames
row.names(tab_prommetaclass)<- row.names(otumetaclass)  #colocando otus en rownames tabla de otus

## Cargando objeto de phyloseq
OTUMETAGclass = otu_table(tab_prommetaclass, taxa_are_rows = TRUE)
TAXMETAGclass = tax_table(otumetaclass); row.names(TAXMETAGclass) <- row.names(otumetaclass); colnames(TAXMETAGclass) <- c("Class")
obphyloseqMETAGclass <- phyloseq(OTUMETAGclass, TAXMETAGclass)

METAGclass <- obphyloseqMETAGclass  %>% 
  tax_glom(taxrank = "Class")  %>%                  # agglomerate at phylum level
  transform_sample_counts(function(x) {100 *x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.5) %>%                         # Filter out low abundance taxa
  arrange(Class)
View(METAGclass)

write.table(METAGclass, file="RelativeAbundanceClassMETAGENOMA.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE) #exportando datos desde R a ubicaci??n local

plotbarsMETAGclass <- ggplot(METAGclass, aes(x = Sample, y = Abundance, fill = Class, order = sort(METAGclass$Class, decreasing = TRUE))) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = phylum_colors) +
  theme(axis.text=element_text(size=15, colour = "black", family = "Times"),           #tama??o del texto de ejes
        axis.title=element_text(size=15, family = "Times"),
        legend.text=element_text(size=15, face = "italic", family = "Times")) +        #tama??o t??tulo eje x
  scale_y_continuous(labels = percent_format()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("\nRelative Abundance") + 
  xlab("\nSites")

plotbarsMETAGclass





###### amp_rarecurve ######

#' Calculate rarefaction curve for each sample.
#'
#' Calculate rarefaction curve for each sample using the vegan rarecurve function directly from a phyloseq object.
#'
#' @usage amp_rarecurve(data)
#'
#' @param data (required) A phyloseq object.
#' @param step Step size for sample sizes in rarefaction curves (default: 100).
#' @param ylim vector of y-axis limits.
#' @param xlim vector of x-axis limits.
#' @param label Label rarefaction curves (default: F).
#' @param color Color lines by metadata.
#' @param color.vector Vector with colors e.g. c("red","white") (default: NULL).
#' @param legend Add a legend to the plot if color is used (default: T).
#' @param legend.position Position of the legend (default: "topleft").
#' 
#' @export
#' @import phyloseq
#' @import vegan
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_rarecurve <- function(data, step = 100, ylim = NULL, xlim = NULL, label = T, color = NULL, legend = T, color.vector = NULL, legend.position = "topleft"){
  
  abund = otu_table(data)@.Data %>% as.data.frame()
  
  if (!is.null(color)) {
    gg_color_hue <- function(n) {
      hues = seq(20, 375, length=n+1)
      hcl(h=hues, l=65, c=100)[1:n]
    }
    group_vector<-sample_data(data)[,color]@.Data %>% as.data.frame()
    names(group_vector)<-"color_variable"
    group_vector<-as.character(group_vector)
    groups<-unique(group_vector)
    n = length(groups)
    cols = gg_color_hue(n)
    if (!is.null(color.vector)){ cols <- color.vector}
    
    col_vector<-rep(brewer.pal(7, "Accent"),length(group_vector))
    for (i in 1:length(group_vector)){
      col_vector[i]<-cols[match(group_vector[i],groups)]
    }
  } else {
    col_vector = brewer.pal(7, "Accent")
    col_vector = sort(met$groups)
  }
  
  if (is.null(ylim) & is.null(xlim)){
    rarecurve(t(abund), step = step, label = label, col = col_vector)
  }
  if (!is.null(ylim) & !is.null(xlim)){
    rarecurve(t(abund), step = step, ylim = ylim, xlim = xlim, label = label, col = col_vector)
  }
  if (!is.null(ylim) & is.null(xlim)){
    rarecurve(t(abund), step = step, ylim = ylim, label = label, col = col_vector)
  }
  if (is.null(ylim) & !is.null(xlim)){
    rarecurve(t(abund), step = step, xlim = xlim, label = label, col = col_vector)
  }
  
  if (!is.null(color) & legend == T){
    legend(legend.position,legend = met$ID_ref,fill = cols, bty = "n")
  }
}