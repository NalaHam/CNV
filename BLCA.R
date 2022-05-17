library(TCGAbiolinks)
library(dplyr)

#BLCA (has normal 1-4 cancer levels)
clin <- GDCquery_clinic("TCGA-BLCA", type = "clinical", save.csv = TRUE)
sub_clin <- subset(clin, select = c("submitter_id", "ajcc_pathologic_stage", 
                                    "gender", "age_at_index", "race", "ethnicity", 
                                    "cigarettes_per_day", "years_smoked","alcohol_history", 
                                    "alcohol_intensity"))

#Somatic Mutations:
query1 <- GDCquery( project = "TCGA-BLCA", 
                    data.category = "Simple Nucleotide Variation", 
                    data.type = "Masked Somatic Mutation", legacy=F)

GDCdownload(query1)

muts <- GDCprepare(query1)

mut_results <- getResults(query1) #can skip

muts <- subset(muts, select = c( "Hugo_Symbol", "Chromosome", "Start_Position", 
                                 "End_Position", "Variant_Classification",
                                 "Variant_Type", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
                                 "Mutation_Status", "HGVSp"))


muts$submitter_id <- substr(muts$Tumor_Sample_Barcode, 1, 12)

mut_data <- merge(muts, sub_clin, by = "submitter_id", all = T)

names(mut_data)[1] <- 'sample_id'
names(mut_data)[2] <- 'GeneSymbol'
names(mut_data)[12] <- 'cancer_level'

write.csv(mut_data,"BLCA_Mut.csv")


#separate into stages and ages
mut_data_1 <- filter(mut_data, cancer_level == "Stage I" | cancer_level == "Stage II", 
                     age_at_index <= 55, 
                     na.rm = TRUE)


#Mutations and Tumor Suppressors compared:
tsg <- semi_join(mut_data_1, Human_TSGs, by = "GeneSymbol") #make new df where strings in column GeneSymbol match in df mut_data and df Human_TSGs

#Get information for just X chromosome mutations:
xchr_tsg <- filter(tsg, tsg$Chromosome == "chrX") #make new df for rows in Chromosome column = chrX 

unq <- unique(xchr_tsg[c("GeneSymbol")])

dif_tsg <- cbind(dif_tsg$GeneSymbol, unq$GeneSymbol)

#sums of male mutations
sum(xchr_tsg$GeneSymbol == "KDM6A" & xchr_tsg$gender == "male") #7
sum(xchr_tsg$GeneSymbol == "XIST" & xchr_tsg$gender == "male") #1
sum(xchr_tsg$GeneSymbol == "RBBP7" & xchr_tsg$gender == "male") #0

unq$male <- c(7,1,0)

#sums of female mutations
sum(xchr_tsg$GeneSymbol == "KDM6A" & xchr_tsg$gender == "female") #1
sum(xchr_tsg$GeneSymbol == "XIST" & xchr_tsg$gender == "female") #0
sum(xchr_tsg$GeneSymbol == "RBBP7" & xchr_tsg$gender == "female") #1

unq$female <- c(1,0,1)

#find difference btw males and females for each tumor suppressor
blca_tsg$blca_bias <- unq$male - unq$female
write.csv(blca_tsg,"BLCA_ratio.csv")




