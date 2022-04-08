library(TCGAbiolinks)
library(dplyr)

#KIRC (has normal 1-4 cancer levels)
clin <- GDCquery_clinic("TCGA-KIRC", type = "clinical", save.csv = TRUE)
sub_clin <- subset(clin, select = c("submitter_id", "ajcc_pathologic_stage", 
                                    "gender", "age_at_index"))

#Somatic Mutations:
query1 <- GDCquery( project = "TCGA-KIRC", 
                    data.category = "Simple Nucleotide Variation", 
                    data.type = "Masked Somatic Mutation", legacy=F)

GDCdownload(query1)

muts <- GDCprepare(query1)

mut_results <- getResults(query1)

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

#separate into stages and ages
mut_data_1 <- filter(mut_data, cancer_level == "Stage I" | cancer_level == "Stage II", 
                     age_at_index <= 55, 
                     na.rm = TRUE)


#Mutations and Tumor Suppressors compared:
tsg <- semi_join(mut_data_1, Human_TSGs, by = "GeneSymbol") #make new df where strings in column GeneSymbol match in df mut_data and df Human_TSGs

#Get information for just X chromosome mutations:
xchr_tsg <- filter(tsg, tsg$Chromosome == "chrX") #make new df for rows in Chromosome column = chrX 

dif_tsg <- unique(xchr_tsg[c("GeneSymbol")])

#sums of male mutations
sum(xchr_tsg$GeneSymbol == "BCORL1" & xchr_tsg$gender == "male") #2
sum(xchr_tsg$GeneSymbol == "DMD" & xchr_tsg$gender == "male") #0
sum(xchr_tsg$GeneSymbol == "XIST" & xchr_tsg$gender == "male") #1

dif_tsg$male <- c(2,0,1)

#sums of female mutations
sum(xchr_tsg$GeneSymbol == "BCORL1" & xchr_tsg$gender == "female") #3
sum(xchr_tsg$GeneSymbol == "DMD" & xchr_tsg$gender == "female") #1
sum(xchr_tsg$GeneSymbol == "XIST" & xchr_tsg$gender == "female") #1

dif_tsg$female <- c(3,1,1)

#find difference btw males and females for each tumor suppressor
dif_tsg$kirc_bias <- dif_tsg$male - dif_tsg$female







