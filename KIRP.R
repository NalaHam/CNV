#KIRP
clin <- GDCquery_clinic("TCGA-KIRP", type = "clinical", save.csv = TRUE)
sub_clin <- subset(clin, select = c("submitter_id", "ajcc_pathologic_stage", 
                                    "gender", "age_at_index", "race", "ethnicity", 
                                    "cigarettes_per_day", "years_smoked","alcohol_history", 
                                    "alcohol_intensity"))

#Somatic Mutations:
query1 <- GDCquery( project = "TCGA-KIRP", 
                    data.category = "Simple Nucleotide Variation", 
                    data.type = "Masked Somatic Mutation", legacy=F)

GDCdownload(query1)

muts <- GDCprepare(query1)

muts <- subset(muts, select = c( "Hugo_Symbol", "Chromosome", "Start_Position", 
                                 "End_Position", "Variant_Classification",
                                 "Variant_Type", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
                                 "Mutation_Status", "HGVSp"))


muts$submitter_id <- substr(muts$Tumor_Sample_Barcode, 1, 12)

mut_data <- merge(muts, sub_clin, by = "submitter_id", all = T)

names(mut_data)[1] <- 'sample_id'
names(mut_data)[2] <- 'GeneSymbol'
names(mut_data)[12] <- 'cancer_level'

write.csv(mut_data,"KIRP_Mut.csv")