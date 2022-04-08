library(TCGAbiolinks)
library(dplyr)

#GBM (has 0-4 cancer levels, BUT NO DATA ON THIS!!! just NA)
clin <- GDCquery_clinic("TCGA-GBM", type = "clinical", save.csv = TRUE)
sub_clin <- subset(clin, select = c("submitter_id", "ajcc_pathologic_stage", 
                                    "gender", "age_at_index"))

#Somatic Mutations:
query1 <- GDCquery( project = "TCGA-GBM", 
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

write.csv(mut_data,"GBM_Mut.csv")


#separate into stages and ages
mut_data_1 <- filter(mut_data, 
                     age_at_index <= 55, 
                     na.rm = TRUE)


#Mutations and Tumor Suppressors compared:
tsg <- semi_join(mut_data_1, Human_TSGs, by = "GeneSymbol") #make new df where strings in column GeneSymbol match in df mut_data and df Human_TSGs

#Get information for just X chromosome mutations:
xchr_tsg <- filter(tsg, tsg$Chromosome == "chrX") #make new df for rows in Chromosome column = chrX 

gbm_tsg <- unique(xchr_tsg[c("GeneSymbol")])

#sums of male mutations
gbm_tsg$male <- c(sum(xchr_tsg$GeneSymbol == "DMD" & xchr_tsg$gender == "male"), #2
sum(xchr_tsg$GeneSymbol == "XIST" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "BCORL1" & xchr_tsg$gender == "male"), #2
sum(xchr_tsg$GeneSymbol == "MIR509-3" & xchr_tsg$gender == "male"), #0
sum(xchr_tsg$GeneSymbol == "RBBP7" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "SRPX" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "FOXP3" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "AMER1" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "BTK" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "BEX2" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "TREX2" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "FLNA" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "RPL10" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "PHF6" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "RBMX" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "FHL1" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "DUSP9" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "KDM6A" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "DDX3X" & xchr_tsg$gender == "male"), #1
sum(xchr_tsg$GeneSymbol == "MIR424" & xchr_tsg$gender == "male")) #1

#sums of female mutations
gbm_tsg$female <- c(sum(xchr_tsg$GeneSymbol == "DMD" & xchr_tsg$gender == "female"), #2
                  sum(xchr_tsg$GeneSymbol == "XIST" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "BCORL1" & xchr_tsg$gender == "female"), #2
                  sum(xchr_tsg$GeneSymbol == "MIR509-3" & xchr_tsg$gender == "female"), #0
                  sum(xchr_tsg$GeneSymbol == "RBBP7" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "SRPX" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "FOXP3" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "AMER1" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "BTK" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "BEX2" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "TREX2" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "FLNA" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "RPL10" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "PHF6" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "RBMX" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "FHL1" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "DUSP9" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "KDM6A" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "DDX3X" & xchr_tsg$gender == "female"), #1
                  sum(xchr_tsg$GeneSymbol == "MIR424" & xchr_tsg$gender == "female")) #1


#find difference btw males and females for each tumor suppressor
gbm_tsg$gbm_bias <- gbm_tsg$male - gbm_tsg$female






