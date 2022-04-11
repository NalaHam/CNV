library(TCGAbiolinks)
library(dplyr)

#LUAD
clin <- GDCquery_clinic("TCGA-LUAD", type = "clinical", save.csv = TRUE)
sub_clin <- subset(clin, select = c("submitter_id", "ajcc_pathologic_stage", 
                                    "gender", "age_at_index"))

#Somatic Mutations:
query1 <- GDCquery( project = "TCGA-LUAD", 
                    data.category = "Simple Nucleotide Variation", 
                    data.type = "Masked Somatic Mutation", legacy=F)

GDCdownload(query1)

muts <- GDCprepare(query1)

muts_sub <- muts[muts$Variant_Classification != "Silent", ] #remove rows that are silent mutations

muts_sub <- subset(muts_sub, select = c( "Hugo_Symbol", "Chromosome", "Start_Position", 
                                         "End_Position", "Variant_Classification",
                                         "Variant_Type", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
                                         "Mutation_Status", "HGVSp"))


submitter_id <- substr(muts_sub$Tumor_Sample_Barcode, 1, 12)

muts_sub$submitter_id <- submitter_id

mut_data <- merge(muts_sub, sub_clin, by = "submitter_id", all = T)

names(mut_data)[1] <- 'sample_id'
names(mut_data)[2] <- 'GeneSymbol'
names(mut_data)[12] <- 'cancer_level'

write.csv(mut_data,"LUAD_tsg.csv")


#separate into stages and ages
unique(clin[c("ajcc_pathologic_stage")])
mut_data_1 <- filter(mut_data, cancer_level == "Stage I"| cancer_level == "Stage IA"
                     | cancer_level == "Stage IB" | cancer_level == "Stage II" 
                     | cancer_level == "Stage II"| cancer_level == "Stage IIA"
                     | cancer_level == "Stage IIB",
                     age_at_index <= 55, 
                     na.rm = TRUE)

#Mutations and Tumor Suppressors compared:
tsg <- semi_join(mut_data_1, Human_TSGs, by = "GeneSymbol") #make new df where strings in column GeneSymbol match in df mut_data and df Human_TSGs

#Get information for just X chromosome mutations:
xchr_tsg <- filter(tsg, tsg$Chromosome == "chrX") #make new df for rows in Chromosome column = chrX 

luad_tsg <- unique(xchr_tsg[c("GeneSymbol")])

#sums of male mutations
luad_tsg$male <- c(sum(xchr_tsg$GeneSymbol == "DMD" & xchr_tsg$gender == "male"),
                   sum(xchr_tsg$GeneSymbol == "ZNF185" & xchr_tsg$gender == "male"),
                   sum(xchr_tsg$GeneSymbol == "FHL1" & xchr_tsg$gender == "male"),
                   sum(xchr_tsg$GeneSymbol == "AMER1" & xchr_tsg$gender == "male"),
                   sum(xchr_tsg$GeneSymbol == "RPS6KA6" & xchr_tsg$gender == "male"),
                   sum(xchr_tsg$GeneSymbol == "BCORL1" & xchr_tsg$gender == "male"),
                   sum(xchr_tsg$GeneSymbol == "FLNA" & xchr_tsg$gender == "male"), 
                   sum(xchr_tsg$GeneSymbol == "XIST" & xchr_tsg$gender == "male"), 
                   sum(xchr_tsg$GeneSymbol == "FOXP3" & xchr_tsg$gender == "male"), 
                   sum(xchr_tsg$GeneSymbol == "BTK" & xchr_tsg$gender == "male"),
                   sum(xchr_tsg$GeneSymbol == "TREX2" & xchr_tsg$gender == "male"),
                   sum(xchr_tsg$GeneSymbol == "EDA2R" & xchr_tsg$gender == "male"),
                   sum(xchr_tsg$GeneSymbol == "KDM6A" & xchr_tsg$gender == "male"),
                   sum(xchr_tsg$GeneSymbol == "DDX3X" & xchr_tsg$gender == "male"))

#sums of female mutations
luad_tsg$female <- c(sum(xchr_tsg$GeneSymbol == "DMD" & xchr_tsg$gender == "female"),
                     sum(xchr_tsg$GeneSymbol == "ZNF185" & xchr_tsg$gender == "female"),
                     sum(xchr_tsg$GeneSymbol == "FHL1" & xchr_tsg$gender == "female"),
                     sum(xchr_tsg$GeneSymbol == "AMER1" & xchr_tsg$gender == "female"),
                     sum(xchr_tsg$GeneSymbol == "RPS6KA6" & xchr_tsg$gender == "female"),
                     sum(xchr_tsg$GeneSymbol == "BCORL1" & xchr_tsg$gender == "female"),
                     sum(xchr_tsg$GeneSymbol == "FLNA" & xchr_tsg$gender == "female"), 
                     sum(xchr_tsg$GeneSymbol == "XIST" & xchr_tsg$gender == "female"), 
                     sum(xchr_tsg$GeneSymbol == "FOXP3" & xchr_tsg$gender == "female"), 
                     sum(xchr_tsg$GeneSymbol == "BTK" & xchr_tsg$gender == "female"),
                     sum(xchr_tsg$GeneSymbol == "TREX2" & xchr_tsg$gender == "female"),
                     sum(xchr_tsg$GeneSymbol == "EDA2R" & xchr_tsg$gender == "female"),
                     sum(xchr_tsg$GeneSymbol == "KDM6A" & xchr_tsg$gender == "female"),
                     sum(xchr_tsg$GeneSymbol == "DDX3X" & xchr_tsg$gender == "female"))

#find difference btw males and females for each tumor suppressor
luad_tsg$lihc_bias <- luad_tsg$male - luad_tsg$female



                                       
                                       