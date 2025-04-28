#marine diatom ISME revisions and analyses

#load packages and data

library(microbiome)
library(vegetarian) 
library(phyloseq); packageVersion("phyloseq") #1.42.0
library(ggplot2); packageVersion("ggplot2") #3.5.1
library(gdata)
library(ecodist)
library(vegan)
library("car") 
library(dplyr)
library(biomformat)
library("ape") #phylogenetic tools packages
library(phytools) #phylogenetic tools packages
library(castor)
library(doParallel) #used for UniFrac but meh not super necessary for small data sets.
library(lme4)

#set seed
set.seed(45)

setwd("/Users/lab/Documents/Cyrus_Daruwala/Microplastics2/Final_Qiime2_outputs")
#to load in the global environment and skip lines 26-248 run the line below.
#load(file="CaddiesSeqAnalysis09262024.rda")

#Read in biom table from working directory
setwd("/Users/lab/Documents/Cyrus_Daruwala/Microplastics2/Final_Qiime2_outputs/mp2_paired_dadatableFINAL")
ASV_reads<-read_biom("feature-table.biom")
ASV_table<-as.data.frame(as.matrix(biom_data(ASV_reads)))
ASV_table[1:10,1:10]
otu_tab<-t(ASV_table)
dim(otu_tab) #215 x 4310
otu_tab[1:10, 1:10]
rownames(otu_tab) #ordered as JRD_CED, 10, 11, 12, etc.. Sample 132 missing due to low read depth

#Read in metadata file
setwd("/Users/lab/Documents/Cyrus_Daruwala/Microplastics2/metadata")
meta_data<-read.csv("Metadata_Final.csv",header=TRUE) #40 obs x 8 vars
head(meta_data)
str(meta_data)
colnames(meta_data)
meta_data<-meta_data[order(meta_data$Initials_SampleNumber),] 
dim(meta_data) #40 x 8
meta_data$Initials_SampleNumber #order matches line 31. 
meta_data$Initials_SampleNumber ==rownames(otu_tab) #bien

#Read in taxonomic information. #First go into Excel, open tsv, and resave as csv.
setwd("/Users/lab/Documents/Cyrus_Daruwala/Microplastics2/Final_Qiime2_outputs/mp2_paired_taxaFINAL2")
bac.taxa<-read.csv(file="MP2_Paired_taxaFINAL3.csv")
str(bac.taxa)
rownames(bac.taxa) #whole numbers
colnames(bac.taxa) # "Feature.ID" "Domain"     "Phylum"     "Class"      "Order"      "Family"     "Genus"      "species"
rownames(bac.taxa)<-bac.taxa[,1] #assign the feature IDs to the row names
bac.taxa<-bac.taxa[,2:8] #get rid of the first column. 
colnames(bac.taxa) # "Domain"  "Phylum"  "Class"   "Order"   "Family"  "Genus"   "species"

#Reading in phylogenetic tree
setwd("/Users/lab/Documents/Cyrus_Daruwala/Microplastics2/Final_Qiime2_outputs/mp2-rooted-tree-paired-seqsFINAL")
dated.16Stree<-read.tree(file="tree.nwk") #reading in tree, uses ape package
is.rooted(dated.16Stree) #TRUE
sample_names(dated.16Stree) #NULL
dated.16Stree$tip.label #for example "41e44b590a0d2c9f9c4e3ea3993981ac" 

#Creating phyloseq object
str(bac.taxa) #data.frame 3814 x 7 
bac.taxa.m<-as.matrix(bac.taxa) #VERY NECESSARY TO DO, DON'T SKIP. 
str(bac.taxa.m)
colnames(bac.taxa.m)
rownames(bac.taxa.m)
str(otu_tab) #37 x 3814
df.bacterial.OTU<-as.data.frame(otu_tab)
str(meta_data) #data.frame 37 x 8

#Matching row names
rownames(df.bacterial.OTU)<-as.character(meta_data[,1])
colnames(df.bacterial.OTU) #accession numbers
rownames(meta_data)<-as.character(meta_data[,1]) #sample names
rownames(bac.taxa.m)<-as.character(rownames(bac.taxa.m)) #these the accession numbers
samp.names<-as.character(meta_data[,1]) 

#Matching sample names (originally marked NULL)
sample_names(df.bacterial.OTU)<-samp.names
sample_names(meta_data)<-samp.names
sample_names(bac.taxa.m)<-samp.names
sample_names(dated.16Stree)<-samp.names

#Matching taxa names (originally marked NULL)
taxa_names(df.bacterial.OTU)<-colnames(df.bacterial.OTU)
taxa_names(dated.16Stree)<-colnames(df.bacterial.OTU)
taxa_names(meta_data)<-colnames(df.bacterial.OTU)
taxa_names(bac.taxa.m)<-colnames(df.bacterial.OTU)

#Here is the actual phyloseq object
Bacterial_phylo<-phyloseq(otu_table(df.bacterial.OTU, taxa_are_rows=FALSE), tax_table(bac.taxa.m), sample_data(meta_data), phy_tree(dated.16Stree)) 

Bacterial_phylo@otu_table[1:10,1:10]
dim(Bacterial_phylo@otu_table) #40 x 3814
rowSums(Bacterial_phylo@otu_table) #
sort(colSums(Bacterial_phylo@otu_table),decreasing=FALSE) # a handful of singletons ##singletons: a sequence that cannot be matched to its counter strand 

#Examining taxonomic ranks to examine where chloroplasts and mitochondria are nested within
table(tax_table(Bacterial_phylo)[, "Domain"], exclude = NULL) #25 unassigned
table(tax_table(Bacterial_phylo)[, "Phylum"], exclude = NULL) #237 unassigned, potentially of interest to remove
table(tax_table(Bacterial_phylo)[, "Class"], exclude = NULL) #examine #244 unassigned
table(tax_table(Bacterial_phylo)[, "Order"], exclude = NULL) #317 reads assigned as chloroplast at this taxonomic rank
#15 uncultured and 284 unassigned
table(tax_table(Bacterial_phylo)[, "Family"], exclude = NULL) #297 reads assigned as mitochondria at this taxonomic rank
table(tax_table(Bacterial_phylo)[, "Genus"], exclude = NULL) #25 unknown family, #683 unassigned
table(tax_table(Bacterial_phylo)[, "Species"], exclude = NULL) #1 unidentified marine bacteria, 2392 unassigned

#Replace taxa names to short hand
#taxa_names(Bacterial_phylo) <- paste0("ASV", seq(ntaxa(Bacterial_phylo)))

#remove other stuff
p1<-subset_taxa(Bacterial_phylo,  !Domain %in% "Unassigned") #requires dplyr (%in% syntax)
p2<-subset_taxa(p1,  !Domain %in% "Archaea") #be no more!
p3<-subset_taxa(p2,  !Domain %in% "Eukaryota") #sayonara euks!
p4<-subset_taxa(p3,  !Order %in% "Chloroplast") #Chlorplasts be gone! 
p5<-subset_taxa(p4,  !Family %in% "Mitochondria") #Adios amigos!

# #double check removal
# table(tax_table(p5)[, "Domain"], exclude = NULL)
# table(tax_table(p5)[, "Phylum"], exclude = NULL) #could be beneficial to blastN the unassigned ones
# table(tax_table(p5)[, "Family"], exclude = NULL)
table(tax_table(p5)[, "Order"], exclude = NULL)


#time to remove the blanks and all the ASVs found in them from all of our samples. 
rowSums(Bacterial_phylo@otu_table) #this will show you read depth before taxonomic filtering
rowSums(p5@otu_table) #this will show you read depth after taxonomic filtering 
dim(p5@otu_table) #215 x #### the difference is the number or taxa aka ASVs assigned to said taxa that you dont care about anuymore

# PCRBlank1<-subset_samples(p5, Treatment.Type=="BLANK1")
# PCRBlank1.rm<-prune_taxa(taxa_sums(PCRBlank1) > 0, PCRBlank1)
# bad_taxa1<-colnames(PCRBlank1.rm@otu_table) #length = 9
# 
# pop_taxa = function(physeq, badTaxa){
#   allTaxa = taxa_names(physeq)
#   allTaxa <- allTaxa[!(allTaxa %in% badTaxa)] #requires dplyr (%in% syntax)
#   return(prune_taxa(allTaxa, physeq))
# }
# 
# dim(p5@otu_table) #36 x 3140
# p6 = pop_taxa(p5, bad_taxa1)
# dim(p6@otu_table) #36 x 3131; #9 less
# str(p6@tax_table)
# str(p6@phy_tree) #works!
# rowSums(p6@otu_table)
# 
# PCRBlank2<-subset_samples(p6, Treatment.Type=="BLANK2")
# PCRBlank2.rm<-prune_taxa(taxa_sums(PCRBlank2) > 0, PCRBlank2)
# bad_taxa2<-colnames(PCRBlank2.rm@otu_table) #length = 4
# 
# dim(p6@otu_table) #40 x 3131
# p7 = pop_taxa(p6, bad_taxa2)
# dim(p7@otu_table) #36 x 3127; #4 less
# rowSums(p7@otu_table)
# 
# PCRBlank3<-subset_samples(p6, Treatment.Type=="BLANK3")
# PCRBlank3.rm<-prune_taxa(taxa_sums(PCRBlank3) > 0, PCRBlank3)
# bad_taxa3<-colnames(PCRBlank3.rm@otu_table) #length = 13
# 
# dim(p7@otu_table) #40 x 3131
# p8 = pop_taxa(p7, bad_taxa3)
# dim(p8@otu_table) #37 x 3116; #13 less
# rowSums(p8@otu_table)

#Remove samples with low read depth
#CaddieD_phylo<-subset_samples(p5, Initials_SampleNumber != "JRD_CED132" & Initials_SampleNumber != "JRD_CED140" & Initials_SampleNumber != "JRD_CED165") & Initials_SampleNumber != "JRD_CED166" & Initials_SampleNumber != "JRD_CED168" & Initials_SampleNumber != "JRD_CED199" & Initials_SampleNumber != "JRD_CED6" & Initials_SampleNumber != "JRD_CED94" & Initials_SampleNumber != "JRD_CED26"
CaddieD_phylo <- subset_samples(p5, 
                                Initials_SampleNumber != "JRD_CED132" & 
                                  Initials_SampleNumber != "JRD_CED140" & 
                                  Initials_SampleNumber != "JRD_CED165" & 
                                  Initials_SampleNumber != "JRD_CED166" & 
                                  Initials_SampleNumber != "JRD_CED168" & 
                                  Initials_SampleNumber != "JRD_CED199" & 
                                  Initials_SampleNumber != "JRD_CED6" & 
                                  Initials_SampleNumber != "JRD_CED94" & 
                                  Initials_SampleNumber != "JRD_CED26")

CaddieD_phylo<-subset_taxa(CaddieD_phylo,  !Phylum %in% "Unassigned") #be no more!
CaddieD_phylo.rm<-prune_taxa(taxa_sums(CaddieD_phylo) > 0, CaddieD_phylo) #this is the object that we want to work with


dim(CaddieD_phylo.rm@sam_data) #37 x 8
CaddieD_meta_pruned<-CaddieD_phylo.rm@sam_data

dim(CaddieD_phylo.rm@otu_table) #37 x 2912

#ASV table with out zeros and singletons ## I dont think this part of the code is significant for me since my read depths were good
CaddieD_phylo.rm2<- prune_taxa(taxa_sums(CaddieD_phylo.rm) > 1, CaddieD_phylo.rm) #remove singletons
dim(CaddieD_phylo.rm2@otu_table) #37 x 2895
obj45<-rowSums(CaddieD_phylo.rm2@otu_table)

CaddieD_phylo.rm2.ASVtable<-CaddieD_phylo.rm2@otu_table

#Rarefying
##accounting for uneven sampling efforts and biases during the sampling process such as if one sample gets amplifyed a lot more for some reason compared to other samples
##will rarefy for the lowest number of reads in our case 
min(rowSums(CaddieD_phylo.rm2.ASVtable)) #28484 
mean(rowSums(CaddieD_phylo.rm2.ASVtable)) #56,383.22
median(rowSums(CaddieD_phylo.rm2.ASVtable)) #51,430
sort(x=rowSums(CaddieD_phylo.rm2.ASVtable), decreasing=FALSE)
bac.tab.df<-as.data.frame(CaddieD_phylo.rm2.ASVtable)
rdat<-rrarefy(bac.tab.df,28484) #rarefy! to minimum number of read depths. 

#Standardize abundances into proportions. 
rowSums(rdat) #Check if it worked
std.bac.tab<-decostand(rdat,"total") #replace object to rdat after rarefying. 
std.bac.tab[1:10,1:10] #Looks groovy! 
rowSums(std.bac.tab)

#Smush back together into single phyloseq object; this contains rarefied and standardized data
bacterial.phylo.4analysis<-phyloseq(otu_table(std.bac.tab, taxa_are_rows=FALSE), sample_data(CaddieD_phylo.rm2@sam_data), tax_table(CaddieD_phylo.rm2@tax_table), phy_tree(CaddieD_phylo.rm2@phy_tree))

#Create Quant Jaccard distance matrix
drdat<-vegdist(std.bac.tab,"bray")

# #metadata
# str(bacterial.phylo.4analysis@sam_data)
# bacterial.phylo.4analysis@otu_table[1:10,1:10]
# bacterial.phylo.4analysis@sam_data$Treatment.Type[bacterial.phylo.4analysis@sam_data$Treatment.Type=="treated"]<-"Antibiotic Treated"
# bacterial.phylo.4analysis@sam_data$Treatment.Type[bacterial.phylo.4analysis@sam_data$Treatment.Type=="untreated"]<-"Untreated"
# bacterial.phylo.4analysis@sam_data$Treatment.Type<-factor(bacterial.phylo.4analysis@sam_data$Treatment.Type, levels=c("Antibiotic Treated","Untreated"))
# #bacterial.phylo.4analysis@sam_data$Filter.type<-factor(bacterial.phylo.4analysis@sam_data$Filter.type)


poop<-bacterial.phylo.4analysis@sam_data
str(poop)
colnames(poop)
poop$Fish<-factor(poop$Fish, levels = c("Fish Present","Fish Absent"))
poop$Nutrients<-factor(poop$Nutrients, levels = c("Increased Nutrients","Ambient Nutrients"))
poop$Plastic<-factor(poop$Plastic, levels = c("Elastollan","TPU181","No Plastic"))
poop$Trt_code<-factor(poop$Trt_code)

basic.mod<-dbrda(drdat~poop$Fish+poop$Nutrients+poop$Plastic)
basic.mod2<-dbrda(drdat~poop$Trt_code)

h<-how(nperm=10000)
anova(basic.mod, permutations = h, by="margin")
# Permutation test for dbrda under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 10000
# 
# Model: dbrda(formula = drdat ~ poop$Fish + poop$Nutrients + poop$Plastic)
# Df SumOfSqs      F    Pr(>F)    
# poop$Fish        1    0.618 1.7725    0.0277 *  
#   poop$Nutrients   1    0.323 0.9260    0.5151    
# poop$Plastic     2    2.729 3.9145 9.999e-05 ***
#   Residual       202   70.421                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1                 

anova(basic.mod2, permutations = h, by="margin")
# Permutation test for dbrda under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 10000
# 
# Model: dbrda(formula = drdat ~ poop$Trt_code)
# Df SumOfSqs      F    Pr(>F)    
# poop$Trt_code  11    7.507 1.9986 9.999e-05 ***
#   Residual      195   66.584                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

wdistUnif<-UniFrac(physeq = bacterial.phylo.4analysis, weighted = TRUE, parallel = TRUE, fast = TRUE)
wUnimod<-dbrda(wdistUnif~poop$Fish+poop$Nutrients+poop$Plastic) #this model matches 177 on the right side of the ~
#wUnimod<-dbrda(wdistUnif~Caddie_meta$Treatment.Type) #this model matches 177 on the right side of the ~

h<-how(nperm=10000)
anova(wUnimod,permutations=h,by="margin")
#Model: dbrda(formula = wdistUnif ~ Caddie_meta$Treatment.Type)
#Df SumOfSqs      F    Pr(>F)    
#Caddie_meta$Treatment.Type  1  0.45917 21.872 9.999e-05 ***
#  Residual                   35  0.73477                     

uwdistUnif<-UniFrac(physeq = bacterial.phylo.4analysis, weighted = FALSE, parallel = TRUE, fast = TRUE)
uwUnimod<-dbrda(uwdistUnif~poop$Fish+poop$Nutrients+poop$Plastic)
h<-how(nperm=10000)
anova(uwUnimod,permutations=h,by="margin")



#Unconstrained Ordinations 
bacterial.phylo.4analysis@sam_data$Date[bacterial.phylo.4analysis@sam_data$Date=="5/28/24"]<-"05/28/24"
bacterial.phylo.4analysis@sam_data$Date[bacterial.phylo.4analysis@sam_data$Date=="6/11/24"]<-"06/11/24"
bacterial.phylo.4analysis@sam_data$Date[bacterial.phylo.4analysis@sam_data$Date=="6/25/24"]<-"06/25/24"
bacterial.phylo.4analysis@sam_data$Date[bacterial.phylo.4analysis@sam_data$Date=="7/17/24"]<-"07/17/24"
bacterial.phylo.4analysis@sam_data$Date[bacterial.phylo.4analysis@sam_data$Date=="8/6/24"]<-"08/06/24"
bacterial.phylo.4analysis@sam_data$Date[bacterial.phylo.4analysis@sam_data$Date=="8/27/24"]<-"08/27/24"
bacterial.phylo.4analysis@sam_data$Date<-factor(bacterial.phylo.4analysis@sam_data$Date,levels=c("05/28/24","06/11/24", "06/25/24", "07/17/24","08/06/24", "08/27/24"))
bacterial.phylo.4analysis@sam_data$Plastic<-factor(bacterial.phylo.4analysis@sam_data$Plastic, level=c("Elastollan","TPU181","No Plastic"))
bacterial.phylo.4analysis@sam_data$Nutrients<-factor(bacterial.phylo.4analysis@sam_data$Nutrients, level=c("Increased Nutrients","Ambient Nutrients"))
bacterial.phylo.4analysis@sam_data$Fish<-factor(bacterial.phylo.4analysis@sam_data$Fish,levels=c("Fish Absent","Fish Present"))
bacterial.phylo.4analysis@sam_data$Trt_code<-factor( bacterial.phylo.4analysis@sam_data$Trt_code)

hfa_bac_pcoa <- ordinate(physeq = bacterial.phylo.4analysis, method = "PCoA", distance = "bray")
hfa_jac_pcoa <- ordinate(physeq = bacterial.phylo.4analysis, method = "PCoA", distance = "jaccard") 
wu_hfa_bac_pcoa <- ordinate(physeq = bacterial.phylo.4analysis, method = "PCoA", distance = "wunifrac") 
uwu_hfa_bac_pcoa <- ordinate(physeq = bacterial.phylo.4analysis, method = "PCoA", distance = "unifrac")

#Plot bray curtis pcoa 
bc1<-plot_ordination(physeq = bacterial.phylo.4analysis, ordination = hfa_bac_pcoa, shape= "Plastic", color = "Date", axes = c(1,2)) + 
  scale_color_manual(values = c("#e95e50","#f29222","#fac723","#a1c65d","#0cb2af","#936fac")) +
  scale_shape_manual(values = c(22,15,8)) +
  #scale_x_continuous(limits=c(-0.5,0.5),breaks=c(-0.5,-0.25,0.0,0.25,0.5))+
  #scale_y_continuous(limits=c(-0.55,0.50),breaks=c(-0.5,-0.25,0.0,0.25,0.5))+
  geom_point(aes(color = Date), alpha = 1, size = 6, stroke=1)+
  theme_test()+
  #labs(color="Treatment Type")+ labs(shape = "Treatment Type") + labs(fill = "Treatment Type") +
  #labs(x="Axis 1 [63.1%]",y="Axis 2 [9.9%]")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right")
bc1$layers<-bc1$layers[-1]
bc1

bc1<-plot_ordination(physeq = bacterial.phylo.4analysis, ordination = hfa_bac_pcoa, shape= "Fish", color = "Date", axes = c(1,2)) + 
  scale_color_manual(values = c("#e95e50","#f29222","#fac723","#a1c65d","#0cb2af","#936fac")) +
  scale_shape_manual(values = c(22,15)) +
  #scale_x_continuous(limits=c(-0.5,0.5),breaks=c(-0.5,-0.25,0.0,0.25,0.5))+
  #scale_y_continuous(limits=c(-0.55,0.50),breaks=c(-0.5,-0.25,0.0,0.25,0.5))+
  geom_point(aes(color = Date), alpha = 1, size = 6, stroke=1)+
  theme_test()+
  #labs(color="Treatment Type")+ labs(shape = "Treatment Type") + labs(fill = "Treatment Type") +
  #labs(x="Axis 1 [63.1%]",y="Axis 2 [9.9%]")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right")
bc1$layers<-bc1$layers[-1]
bc1

bc1c<-plot_ordination(physeq = bacterial.phylo.4analysis, ordination = hfa_bac_pcoa, shape= "Fish", color = "Plastic", axes = c(1,2)) + 
  scale_color_manual(values = c("#e95e50","#f29222","#fac723")) +
  scale_shape_manual(values = c(22,15)) +
  #scale_x_continuous(limits=c(-0.5,0.5),breaks=c(-0.5,-0.25,0.0,0.25,0.5))+
  #scale_y_continuous(limits=c(-0.55,0.50),breaks=c(-0.5,-0.25,0.0,0.25,0.5))+
  geom_point(aes(color = Plastic), alpha = 1, size = 6, stroke=1)+
  theme_test()+
  #labs(color="Treatment Type")+ labs(shape = "Treatment Type") + labs(fill = "Treatment Type") +
  #labs(x="Axis 1 [63.1%]",y="Axis 2 [9.9%]")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right")
bc1c$layers<-bc1c$layers[-1]
bc1c

#weighted unifrac
bc3<-plot_ordination(physeq = bacterial.phylo.4analysis, ordination = wu_hfa_bac_pcoa, shape = "Treatment.Type", color = "Treatment.Type", axes = c(1,2)) + 
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_color_manual(values =c("black","black"), guide = "none") +
  scale_shape_manual(values = c(22,21)) +
  scale_x_continuous(limits=c(-0.4,0.4),breaks=c(-0.4,-0.2,0.0,0.2,0.4))+
  scale_y_continuous(limits=c(-0.3,0.3),breaks=c(-0.3,-0.15,0.0,0.15,0.3))+
  geom_point(aes(fill = Treatment.Type), alpha = 1, size = 6, stroke=1)+
  theme_test()+
  labs(color="Treatment Type")+ labs(shape = "Treatment Type") + labs(fill = "Treatment Type") +
  labs(x="Axis 1 [48.8%]",y="Axis 2 [28.8%]")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 
bc3$layers<-bc3$layers[-1]


#unweighted unifrac
bc4<-plot_ordination(physeq = bacterial.phylo.4analysis, ordination = uwu_hfa_bac_pcoa, shape = "Treatment.Type", color = "Treatment.Type", axes = c(1,2)) + 
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_color_manual(values =c("black","black"), guide = "none") +
  scale_shape_manual(values = c(22,21)) +
  scale_x_continuous(limits=c(-0.3,0.3),breaks=c(-0.3,-0.15,0.0,0.15,0.3))+
  scale_y_continuous(limits=c(-0.3,0.3),breaks=c(-0.3,-0.15,0.0,0.15,0.3))+
  geom_point(aes(fill = Treatment.Type), alpha = 1, size = 6, stroke=1) +
  theme_test()+
  labs(color="Treatment Type")+ labs(shape = "Treatment Type") + labs(fill = "Treatment Type") +
  labs(x="Axis 1 [24.1%]",y="Axis 2 [7.6%]")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 
bc4$layers<-bc4$layers[-1]


#plot unweighted unifrac pcoa
# bc3<-plot_ordination(physeq = bacterial.phylo.4analysis, ordination = uwu_hfa_bac_pcoa, color = "Water_Source", shape = "Timepoint", axes = c(1,2)) + 
#   scale_fill_manual(values = c("#060282","#5aafed","#ffc001","#ff0000")) +
#   scale_color_manual(values =c("black","black","black","black")) +
#   scale_shape_manual(values = c(21,22,23,24,25)) +
#   scale_x_continuous(limits=c(-0.3,0.3),breaks=c(-0.3,-0.15,0.0,0.15,0.3))+
#   #scale_y_continuous(limits=c(-0.3,0.4),breaks=c(-0.3,-0.1,0.0,0.1,0.4))+
#   geom_point(aes(fill = Water_Source), alpha = 1, size = 4, stroke=1)+
#   theme_test()+
#   labs(shape="Time Point",fill="Origin")+
#   labs(x="Axis 1 [16.3%]",y="Axis 2 [14.9%]")+
#   theme(text = element_text(size = 15))+
#   theme(legend.position = "none")
# bc3$layers<-bc3$layers[-1]

#save(list=ls(),file="AntiB_Caddie_16SAnalyses.rda")


### PERMANOVA TUTORIAL ###
# Calculate bray curtis distance matrix
# caddie_bray <- phyloseq::distance(bacterial.phylo.4analysis, method = "bray")
# 
# # make a data frame from the sample_data
# caddie_sampledf <- data.frame(sample_data(Caddie_meta))
# str(caddie_sampledf)
# adonis2(caddie_bray ~ Treatment.Type, data = caddie_sampledf)

#Permutation test for adonis under reduced model
#adonis2(formula = caddie_bray ~ Treatment.Type, data = caddie_sampledf)
#               Df SumOfSqs      R2      F Pr(>F)    
#Treatment.Type  1   4.3826 0.61914 56.897  0.001 ***
#Residual       35   2.6959 0.38086                  
#Total          36   7.0785 1.00000  

# Homogeneity of dispersion test
dist.hill<-function(dat,q=1){
  dmat<-array(dim=c(dim(dat)[1],dim(dat)[1]))
  for(i in 1:dim(dmat)[1]){
    for(j in 1:dim(dmat)[1]){
      dmat[i,j]<-d(rbind(dat[i,],dat[j,]),lev="beta",q=q)-1
    }
  }
  res<-as.dist(dmat)
  return(res)
}

ddatq0<-dist.hill(bacterial.phylo.4analysis@otu_table,q=0) #pairwise beta-diversity
caddiebeta <- betadisper(ddatq0, caddie_sampledf$Treatment.Type)
mod.beta<-lm(caddiebeta$distances ~ caddiebeta$group)
Anova(mod.beta)

#**go to new script for figure. 

#Anova Table (Type II tests)

#Response: caddiebeta$distances
#Sum Sq Df F value    Pr(>F)    
#caddiebeta$group 0.061230  1  45.439 8.274e-08 ***
#  Residuals        0.047163 35                      
---

#the adonis test is significant, meaning that we can reject the null hypothesis 
#when treated and untreated are compared they are significantly different 
#furthermore, the betadisper test reveals that the results are also significant
#we can say that the groups have different dispersions 


### ALPHA DIVERSITY TUTORIAL ###
#estimating alpha diversity can be problematic
#subsample the libraries with replacement while standardizing sampling effort 

#Caddie_alpha <- CaddieD_phylo.rm2@otu_table

#min_lib <- min(sample_sums(Caddie_alpha))

# Initialize matrices to store richness and evenness estimates
#nsamp = nsamples(Caddie_meta)
#trials = 100

#richness <- matrix(nrow = nsamp, ncol = trials)
#row.names(richness) <- sample_names(Caddie_meta)

#evenness <- matrix(nrow = nsamp, ncol = trials)
#row.names(evenness) <- sample_names(Caddie_meta)

# It is always important to set a seed when you subsample so your result is replicable 
#set.seed(3)

#for (i in 1:100) {
  # Subsample
 # caddie_rar <- rarefy_even_depth(Caddie_alpha, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
 # caddie_rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
 # richness[ ,i] <- caddie_rich
  
  # Calculate evenness
 # caddie_even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
#  evenness[ ,i] <- caddie_even
#}

#Let’s calculate the mean and standard deviation per sample for observed richness and inverse simpson’s index and store those values in a dataframe.

# Create a new dataframe to hold the means and standard deviations of richness estimates
#CaddieID <- row.names(richness)
#mean <- apply(richness, 1, mean)
#measure <- rep("Richness", nsamp)
#caddie_rich_stats <- data.frame(CaddieID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
#CaddieID <- row.names(evenness)
#mean <- apply(evenness, 1, mean)
#sd <- apply(evenness, 1, sd)
#measure <- rep("Inverse Simpson", nsamp)
#caddie_even_stats <- data.frame(CaddieID, mean, sd, measure)

#combine richness and evenness into one dataframe 
#caddie_rich_even <- rbind(caddie_rich_stats, caddie_even_stats)
#add the sample metadata into the dataframe using the merge() command 
#s <- data.frame(sample_data(Caddie_meta))
#alphadiv <- merge(caddie_rich_even, s, by = "CaddieID") 

#plotting alpha diversity measures in a timeseries using facet 
#ggplot(alphadiv, aes(x = Date, y = mean, color = Treatment.Type, group = Treatment.Type, shape = Treatment.Type)) +
#  geom_point(size = 2) + 
#  geom_line(size = 0.8) +
#  facet_wrap(~measure, ncol = 1, scales = "free") +
#  scale_color_manual(values = c("#E96446", "#302F3D", "#87CEFA")) +
#  theme(
#    axis.title.x = element_blank(),
#    axis.title.y = element_blank()
#  )


### TRYING ALPHA DIVERSITY A DIFFERENT WAY ###
alpha0 <- apply(rdat, 1, vegetarian::d, q =0)
alpha1 <- apply(rdat, 1, vegetarian::d, q =1)
alpha2 <- apply(rdat, 1, vegetarian::d, q =2)

#attaching alpha diversity values to the meta data 
caddie_alpha_meta <- cbind(Caddie_meta, alpha0, alpha1, alpha2)
#creating a linear model where the effects of whether or not antibiotics were added
caddie_lm_alpha <- lm(alpha0 ~ Treatment.Type, dat = caddie_alpha_meta)
caddie_lm_alpha1<- lm(alpha1 ~ Treatment.Type, dat = caddie_alpha_meta)
caddie_lm_alpha2<- lm(alpha2 ~ Treatment.Type, dat = caddie_alpha_meta)

#running ANOVA for the model just created 
Anova(caddie_lm_alpha) #richness
Anova(caddie_lm_alpha1) #Exponential Shannon
Anova(caddie_lm_alpha2) #Inverse simpsons index

#creating a boxplot figure 
a0.plot <- ggplot(data = caddie_alpha_meta, aes(x = Treatment.Type, y = alpha0, fill = Treatment.Type)) +
  geom_boxplot(aes(fill = Treatment.Type), alpha = 0.65, outlier.size = 0) +
  geom_point(pch = 21, size=2.5, position = position_jitterdodge())+
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_y_continuous(limits=c(150,350),breaks=c(150,200,250,300,350))+
  theme_classic()+
  labs(y = "Gut Microbiome ASV Richness", x = "Treatment Type") +
  theme(legend.position = "none")

a0.plota<-a0.plot + labs(fill = "Treatment Type") + theme(text= element_text(size=15), legend.text.align = 0) 

a1.plot <- ggplot(data = caddie_alpha_meta, aes(x = Treatment.Type, y = alpha1, fill = Treatment.Type)) +
  geom_boxplot(aes(fill = Treatment.Type), alpha = 0.65, outlier.size = 0) +
  geom_point(pch = 21, size=2.5, position = position_jitterdodge())+
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_y_continuous(limits=c(0,30),breaks=c(0,10,20,30))+
  theme_classic()+
  labs(y = "Gut Microbiome Expontential Shannon's Diversity", x = "Treatment Type") +
  theme(legend.position = "none")

a1.plota<-a1.plot + labs(fill = "Treatment Type") + theme(text= element_text(size=15), legend.text.align = 0) 


a2.plot <- ggplot(data = caddie_alpha_meta, aes(x = Treatment.Type, y = alpha2, fill = Treatment.Type)) +
  geom_boxplot(aes(fill = Treatment.Type), alpha = 0.65, outlier.size = 0) +
  geom_point(pch = 21, size=2.5, position = position_jitterdodge())+
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_y_continuous(limits=c(0,15),breaks=c(0,5,10,15))+
  theme_classic()+
  labs(y = "Gut Microbiome Inverse Simpson's Diversity", x = "Treatment Type") +
  theme(legend.position = "right")

a2.plota<-a2.plot + labs(fill = "Treatment Type") + theme(text= element_text(size=15), legend.text.align = 0) 

library(cowplot)
?plot_grid
#save(list=ls(),file="CaddiesSeqAnalysis09262024.rda")
plot_grid(a0.plota,a1.plota,a2.plota,nrow=1,ncol=3,rel_widths=c(1,1,1.5))

#CaddieD_phylo.rm2
depth<-rowSums(CaddieD_phylo.rm2@otu_table)
mmmmme<-CaddieD_phylo.rm2@sam_data
klklklk<-cbind(mmmmme,depth)
str(klklklk)
klklklk$Treatment.Type[klklklk$Treatment.Type=="treated"]<-"Antibiotic Treated"
klklklk$Treatment.Type[klklklk$Treatment.Type=="untreated"]<-"Untreated"
klklklk$Treatment.Type<-factor(klklklk$Treatment.Type, levels=c("Antibiotic Treated","Untreated"))

rd.plot <- ggplot(data = klklklk, aes(x = Treatment.Type, y = depth, fill = Treatment.Type)) +
  geom_boxplot(aes(fill = Treatment.Type), alpha = 0.65, outlier.size = 0) +
  geom_point(pch = 21, size=2.5, position = position_jitterdodge())+
  scale_fill_manual(values = c("#d81159","#b5e2fa")) +
  scale_y_continuous(limits=c(25000,125000),breaks=c(25000,50000,75000,100000,125000))+
  theme_classic()+
  labs(y = "Gut Microbiome Read Depth", x = "Treatment Type") +
  theme(legend.position = "right")

rd.plota<-rd.plot + labs(fill = "Treatment Type") + theme(text= element_text(size=15), legend.text.align = 0) 

klklklk.mod <- lm(depth ~ Treatment.Type, dat = klklklk)
Anova(klklklk.mod)
# Sum Sq Df F value   Pr(>F)    
# Treatment.Type 4.5679e+09  1  14.989 0.000452 ***
#   Residuals      1.0667e+10 35         

klklklk.antiB<-klklklk[which(klklklk$Treatment.Type=="Antibiotic Treated"),]
klklklk.control<-klklklk[which(klklklk$Treatment.Type=="Untreated"),]

mean(klklklk.antiB$depth) #45,468.42
mean(klklklk.control$depth) #67,798.83

dim(klklklk.antiB) #19 x 8
dim(klklklk.control) #18 x 8

sd(klklklk.antiB$depth) #10,939.47
sd(klklklk.control$depth) #22,377.15

10939.47/(sqrt(19)) #2509.687
22377.15/(sqrt(18)) #5274.345

#venn diagram for treatments using rarefied phyloseq object , bacterial.phylo.4analysis
library(ggvenn)
#pull out your two treats
antiB_cadd<-subset_samples(bacterial.phylo.4analysis, Treatment.Type=="Antibiotic Treated")
antiB_cadd.rm<-prune_taxa(taxa_sums(antiB_cadd) > 0, antiB_cadd)

Cont_cadd<-subset_samples(bacterial.phylo.4analysis, Treatment.Type=="Untreated")
Cont_cadd.rm<-prune_taxa(taxa_sums(Cont_cadd) > 0, Cont_cadd)


#pull out just the colnames
lowASV<-colnames(antiB_cadd.rm@otu_table)
medASV<-colnames(Cont_cadd.rm@otu_table)

#make list and plot fig
new_list<-list("Antibiotic Treated" = lowASV, "Untreated" = medASV)
ggvenn(new_list, fill_color = c("#d81159","#b5e2fa"), fill_alpha = 0.8, stroke_size = 0.8, set_name_size = 6, text_size = 6, show_percentage = FALSE)
