# Installation ------------------------------------------------------------
# source("http://bioconductor.org/biocLite.R")
# biocLite("devtools")    # only if devtools not yet installed
# biocLite("COMBINE-lab/wasabi")
#bioclite(')
# biocLite("rhdf5")
# devtools::install_github("pachterlab/sleuth")
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")

# Library -----------------------------------------------------------------
print('Loading Libraries')
library(wasabi)
library('cowplot')
suppressMessages({library("sleuth")})
library('biomaRt')

# Prepare Samples for Sleuth ----------------------------------------------
#setwd('current_directory')

#only need run once
#sample_id<-read.table(file='./salmon/directory.txt')
#sfdirs <- file.path("./salmon", as.character(sample_id$V1))
#prepare_fish_for_sleuth(sfdirs)

# Sleuth ------------------------------------------------------------------
print('Reading In Sample Info')
sample_id <- dir(file.path("../Results/", "salmon")) #sample_ids
s2c <- read.table(file.path("../", "Sample_Sheet.csv"), header = TRUE,sep=',', stringsAsFactors=FALSE) #metadata
colnames(s2c)[1]<-'sample'

# Ensembl Transcript Mapping Grab -----------------------------------------

print("Gene to Transcript Mapping")

# #determine which ensembl database to use
# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# 
# #grab gene/transcript mapping
# gt_mapping<- getBM(attributes = c("ensembl_transcript_id", "transcript_version",
#                                   "ensembl_gene_id", "external_gene_name", "description",
#                                   "transcript_biotype"),mart = ensembl)
# 
# #rename columns and select mapping
# gt_mapping <- dplyr::rename(gt_mapping, transcript_id = ensembl_transcript_id,
#                             gene_id = ensembl_gene_id, gene_name = external_gene_name)
# 
# #load in data
# gt_mapping <- dplyr::rename(gt_mapping, target_id = transcript_id,ens_gene = gene_id, ext_gene = gene_name)
# 
# save(gt_mapping,file="gt_mapping.Rda")

load("gt_mapping.Rda")

# Function ----------------------------------------------------------------
analyze_sleuth_subset_multi_D5<-function(Variable,Evaluation,Separate,analyze,name){
  print('Sample Information')
  #sample & condition (multiple variables for the diff conditions to be tested)
  
  s2c_treatment<- dplyr::select(subset(s2c,Day=='D5'),
                                sample = description, condition=Treatment,Donor,Day)
  
  
  #s2c_treatment<- dplyr::select(subset(s2c,eval(parse(text=paste(Variable,Evaluation,Separate)))),
  #                             sample = description, condition=Treatment,Donor,Day)
 
  #filter for samples quantified
  s2c_treatment <- s2c_treatment[s2c_treatment$sample %in% sample_id,]
  #print('Samples selected & Condition')
  #print(s2c_treatment)
  
  #add path to samples
  #s2c_treatment <- dplyr::mutate(s2c_treatment, path = file.path("./salmon", sample_id[sample_id %in% s2c_treatment$sample] ))
  s2c_treatment <- dplyr::mutate(s2c_treatment, path = file.path("../Results/salmon", sample_id[sample_id %in% s2c_treatment$sample] ))
  
  #sleuth run
  print('Running Sleuth')
  so_treatment <- sleuth_prep(s2c_treatment,target_mapping = gt_mapping,
                              aggregation_column = 'ens_gene',extra_bootstrap_summary = TRUE)
  
  #statistics [TODO FIX]
  print('Running Statistics')
  so_treatment <- sleuth_fit(so_treatment, ~condition, 'full')
  so_treatment <- sleuth_fit(so_treatment, ~1, 'reduced')
  so_treatment <- sleuth_wt(so_treatment, 'conditionOLA', 'full')
  
  #so_treatment <- sleuth_fit(so_treatment, ~condition + Day + Donor, 'full')
  #so_treatment <- sleuth_fit(so_treatment, ~condition + Day, 'condition_day')
  #so_treatment <- sleuth_fit(so_treatment, ~Day + Donor, 'day_donor')
  #so_treatment <- sleuth_lrt(so_treatment, 'reduced', 'full')
  #results_wt <- sleuth_wt(so_treatment, 'conditionOLA', 'full')
  
  #tables
  #print('Producing Tables')
  #so_treatment_sleuth_table <- sleuth_results(results_wt,test = 'conditionOLA',
  #                                            which_model = 'full',test_type = 'wt', show_all = FALSE)
  #so_treatment_significant <- dplyr::filter(so_treatment_sleuth_table, qval <= 0.05)
  
  #save sleuth objects
  print('Saving Sleuth Object')
  output<-paste(name,'.sleuthobj',sep="")
  sleuth_save(so_treatment, file=output)
  
  #save tables
  #print('Saving Tables')
  #output<-paste(name,'_sleuth_results.tsv',sep="")
  #write.table(so_treatment_sleuth_table,output,quote=FALSE,row.names=FALSE)
  #output<-paste(name,'_sleuth_sig_results.tsv',sep="")
  #write.table(so_treatment_significant,output,quote=FALSE,row.names=FALSE)
}
analyze_sleuth_subset_multi_D6<-function(Variable,Evaluation,Separate,analyze,name){
  print('Sample Information')
  #sample & condition (multiple variables for the diff conditions to be tested)
  
  s2c_treatment<- dplyr::select(subset(s2c,Day=='D6'),
                                sample = description, condition=Treatment,Donor,Day)
  
  
  #s2c_treatment<- dplyr::select(subset(s2c,eval(parse(text=paste(Variable,Evaluation,Separate)))),
  #                             sample = description, condition=Treatment,Donor,Day)
  
  #filter for samples quantified
  s2c_treatment <- s2c_treatment[s2c_treatment$sample %in% sample_id,]
  #print('Samples selected & Condition')
  #print(s2c_treatment)
  
  #add path to samples
  #s2c_treatment <- dplyr::mutate(s2c_treatment, path = file.path("./salmon", sample_id[sample_id %in% s2c_treatment$sample] ))
  s2c_treatment <- dplyr::mutate(s2c_treatment, path = file.path("../Results/salmon", sample_id[sample_id %in% s2c_treatment$sample] ))
  
  #sleuth run
  print('Running Sleuth')
  so_treatment <- sleuth_prep(s2c_treatment,target_mapping = gt_mapping,
                              aggregation_column = 'ens_gene',extra_bootstrap_summary = TRUE)
  
  #statistics [TODO FIX]
  print('Running Statistics')
  so_treatment <- sleuth_fit(so_treatment, ~condition, 'full')
  so_treatment <- sleuth_fit(so_treatment, ~1, 'reduced')
  so_treatment <- sleuth_wt(so_treatment, 'conditionOLA', 'full')
  
  #so_treatment <- sleuth_fit(so_treatment, ~condition + Day + Donor, 'full')
  #so_treatment <- sleuth_fit(so_treatment, ~condition + Day, 'condition_day')
  #so_treatment <- sleuth_fit(so_treatment, ~Day + Donor, 'day_donor')
  #so_treatment <- sleuth_lrt(so_treatment, 'reduced', 'full')
  #results_wt <- sleuth_wt(so_treatment, 'conditionOLA', 'full')
  
  #tables
  #print('Producing Tables')
  #so_treatment_sleuth_table <- sleuth_results(results_wt,test = 'conditionOLA',
  #                                            which_model = 'full',test_type = 'wt', show_all = FALSE)
  #so_treatment_significant <- dplyr::filter(so_treatment_sleuth_table, qval <= 0.05)
  
  #save sleuth objects
  print('Saving Sleuth Object')
  output<-paste(name,'.sleuthobj',sep="")
  sleuth_save(so_treatment, file=output)
  
  #save tables
  #print('Saving Tables')
  #output<-paste(name,'_sleuth_results.tsv',sep="")
  #write.table(so_treatment_sleuth_table,output,quote=FALSE,row.names=FALSE)
  #output<-paste(name,'_sleuth_sig_results.tsv',sep="")
  #write.table(so_treatment_significant,output,quote=FALSE,row.names=FALSE)
}
analyze_sleuth_subset_multi_D9<-function(Variable,Evaluation,Separate,analyze,name){
  print('Sample Information')
  #sample & condition (multiple variables for the diff conditions to be tested)
  
  s2c_treatment<- dplyr::select(subset(s2c,Day=='D9'),
                                sample = description, condition=Treatment,Donor,Day)
  
  
  #s2c_treatment<- dplyr::select(subset(s2c,eval(parse(text=paste(Variable,Evaluation,Separate)))),
  #                             sample = description, condition=Treatment,Donor,Day)
  
  #filter for samples quantified
  s2c_treatment <- s2c_treatment[s2c_treatment$sample %in% sample_id,]
  #print('Samples selected & Condition')
  #print(s2c_treatment)
  
  #add path to samples
  #s2c_treatment <- dplyr::mutate(s2c_treatment, path = file.path("./salmon", sample_id[sample_id %in% s2c_treatment$sample] ))
  s2c_treatment <- dplyr::mutate(s2c_treatment, path = file.path("../Results/salmon", sample_id[sample_id %in% s2c_treatment$sample] ))
  
  #sleuth run
  print('Running Sleuth')
  so_treatment <- sleuth_prep(s2c_treatment,target_mapping = gt_mapping,
                              aggregation_column = 'ens_gene',extra_bootstrap_summary = TRUE)
  
  #statistics [TODO FIX]
  print('Running Statistics')
  so_treatment <- sleuth_fit(so_treatment, ~condition, 'full')
  so_treatment <- sleuth_fit(so_treatment, ~1, 'reduced')
  so_treatment <- sleuth_wt(so_treatment, 'conditionOLA', 'full')
  
  #so_treatment <- sleuth_fit(so_treatment, ~condition + Day + Donor, 'full')
  #so_treatment <- sleuth_fit(so_treatment, ~condition + Day, 'condition_day')
  #so_treatment <- sleuth_fit(so_treatment, ~Day + Donor, 'day_donor')
  #so_treatment <- sleuth_lrt(so_treatment, 'reduced', 'full')
  #results_wt <- sleuth_wt(so_treatment, 'conditionOLA', 'full')
  
  #tables
  #print('Producing Tables')
  #so_treatment_sleuth_table <- sleuth_results(results_wt,test = 'conditionOLA',
  #                                            which_model = 'full',test_type = 'wt', show_all = FALSE)
  #so_treatment_significant <- dplyr::filter(so_treatment_sleuth_table, qval <= 0.05)
  
  #save sleuth objects
  print('Saving Sleuth Object')
  output<-paste(name,'.sleuthobj',sep="")
  sleuth_save(so_treatment, file=output)
  
  #save tables
  #print('Saving Tables')
  #output<-paste(name,'_sleuth_results.tsv',sep="")
  #write.table(so_treatment_sleuth_table,output,quote=FALSE,row.names=FALSE)
  #output<-paste(name,'_sleuth_sig_results.tsv',sep="")
  #write.table(so_treatment_significant,output,quote=FALSE,row.names=FALSE)
}
analyze_sleuth_subset_multi_D6_vs_D9<-function(Variable,Evaluation,Separate,analyze,name){
  print('Sample Information')
  #sample & condition (multiple variables for the diff conditions to be tested)
  
  #s2c_treatment<- dplyr::select(subset(s2c,Treatment=='OLA'),
  #                              sample = description, condition=Treatment,Donor,Day)
 
  s2c_treatment<- dplyr::select(subset(s2c_treatment,Day == 'D6' | Day == 'D9'),
                                sample, condition=Treatment,Donor,Day)   
  
  #s2c_treatment<- dplyr::select(subset(s2c,eval(parse(text=paste(Variable,Evaluation,Separate)))),
  #                             sample = description, condition=Treatment,Donor,Day)
  
  #filter for samples quantified
  s2c_treatment <- s2c_treatment[s2c_treatment$sample %in% sample_id,]
  #print('Samples selected & Condition')
  #print(s2c_treatment)
  
  #add path to samples
  #s2c_treatment <- dplyr::mutate(s2c_treatment, path = file.path("./salmon", sample_id[sample_id %in% s2c_treatment$sample] ))
  s2c_treatment <- dplyr::mutate(s2c_treatment, path = file.path("../Results/salmon", sample_id[sample_id %in% s2c_treatment$sample] ))
  
  #sleuth run
  print('Running Sleuth')
  so_treatment <- sleuth_prep(s2c_treatment,target_mapping = gt_mapping,
                              aggregation_column = 'ens_gene',extra_bootstrap_summary = TRUE)
  
  #statistics [TODO FIX]
  print('Running Statistics')
  so_treatment <- sleuth_fit(so_treatment, ~Day + condition, 'full')
  so_treatment <- sleuth_fit(so_treatment, ~1, 'reduced')
  so_treatment <- sleuth_wt(so_treatment, 'DayD9', 'full')
  
  #so_treatment <- sleuth_fit(so_treatment, ~condition + Day + Donor, 'full')
  #so_treatment <- sleuth_fit(so_treatment, ~condition + Day, 'condition_day')
  #so_treatment <- sleuth_fit(so_treatment, ~Day + Donor, 'day_donor')
  #so_treatment <- sleuth_lrt(so_treatment, 'reduced', 'full')
  #results_wt <- sleuth_wt(so_treatment, 'conditionOLA', 'full')
  
  #tables
  #print('Producing Tables')
  #so_treatment_sleuth_table <- sleuth_results(results_wt,test = 'conditionOLA',
  #                                            which_model = 'full',test_type = 'wt', show_all = FALSE)
  #so_treatment_significant <- dplyr::filter(so_treatment_sleuth_table, qval <= 0.05)
  
  #save sleuth objects
  print('Saving Sleuth Object')
  output<-paste(name,'.sleuthobj',sep="")
  sleuth_save(so_treatment, file=output)
  
  #save tables
  #print('Saving Tables')
  #output<-paste(name,'_sleuth_results.tsv',sep="")
  #write.table(so_treatment_sleuth_table,output,quote=FALSE,row.names=FALSE)
  #output<-paste(name,'_sleuth_sig_results.tsv',sep="")
  #write.table(so_treatment_significant,output,quote=FALSE,row.names=FALSE)
}
analyze_sleuth_subset_olaparib<-function(Variable,Evaluation,Separate,analyze,name){
  print('Sample Information')
  #sample & condition (multiple variables for the diff conditions to be tested)
  s2c_treatment<- dplyr::select(s2c,sample = description, condition=Treatment,Donor,Day)
  
  #add path to samples
  #s2c_treatment <- dplyr::mutate(s2c_treatment, path = file.path("./salmon", sample_id[sample_id %in% s2c_treatment$sample] ))
  s2c_treatment <- dplyr::mutate(s2c_treatment, path = file.path("../Results/salmon", sample_id[sample_id %in% s2c_treatment$sample] ))
  
  #sleuth run
  print('Running Sleuth')
  so_treatment <- sleuth_prep(s2c_treatment,target_mapping = gt_mapping,
                              aggregation_column = 'ens_gene',extra_bootstrap_summary = TRUE)
  
  #statistics [TODO FIX]
  print('Running Statistics')
  so_treatment <- sleuth_fit(so_treatment, ~condition + Day , 'full')
  so_treatment <- sleuth_fit(so_treatment, ~1, 'reduced')
  
  #lrt analysis across all conditions == what transcripts change the most across all conditions
  so_treatment <- sleuth_lrt(so_treatment, 'reduced', 'full')
  
  
  
  #results
  #sleuth_table <- sleuth_results(so_treatment, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)
  #sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
  
  #sleuth_live(so_treatment)
  
  #save sleuth objects
  print('Saving Sleuth Object')
  output<-paste(name,'.sleuthobj',sep="")
  sleuth_save(so_treatment, file=output)
  
  #save tables
  #print('Saving Tables')
  #output<-paste(name,'_sleuth_results.tsv',sep="")
  #write.table(so_treatment_sleuth_table,output,quote=FALSE,row.names=FALSE)
  #output<-paste(name,'_sleuth_sig_results.tsv',sep="")
  #write.table(so_treatment_significant,output,quote=FALSE,row.names=FALSE)
}
analyze_sleuth_subset_D5_D9<-function(Variable,Evaluation,Separate,analyze,name){
  print('Sample Information')
  #sample & condition (multiple variables for the diff conditions to be tested)
  s2c_treatment<- dplyr::select(subset(s2c,Day=='D5' | Day=='D9'),
                                sample = description, condition=Treatment,Donor,Day)
  
  #add path to samples
  #s2c_treatment <- dplyr::mutate(s2c_treatment, path = file.path("./salmon", sample_id[sample_id %in% s2c_treatment$sample] ))
  s2c_treatment <- dplyr::mutate(s2c_treatment, path = file.path("../Results/salmon", sample_id[sample_id %in% s2c_treatment$sample] ))
  
  #sleuth run
  print('Running Sleuth')
  so_treatment <- sleuth_prep(s2c_treatment,target_mapping = gt_mapping,
                              aggregation_column = 'ens_gene',extra_bootstrap_summary = TRUE)
  
  #statistics [TODO FIX]
  print('Running Statistics')
  so_treatment <- sleuth_fit(so_treatment, ~condition + Day , 'full')
  so_treatment <- sleuth_fit(so_treatment, ~Day, 'reduced')
  
  #lrt analysis across all conditions == what transcripts change the most across all conditions
  so_treatment <- sleuth_lrt(so_treatment, 'reduced', 'full')
  
  #results
  #sleuth_table <- sleuth_results(so_treatment, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)
  #sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
  
  #sleuth_live(so_treatment)
  
  #save sleuth objects
  print('Saving Sleuth Object')
  output<-paste(name,'.sleuthobj',sep="")
  sleuth_save(so_treatment, file=output)
  
  #save tables
  #print('Saving Tables')
  #output<-paste(name,'_sleuth_results.tsv',sep="")
  #write.table(so_treatment_sleuth_table,output,quote=FALSE,row.names=FALSE)
  #output<-paste(name,'_sleuth_sig_results.tsv',sep="")
  #write.table(so_treatment_significant,output,quote=FALSE,row.names=FALSE)
}

# Analysis ----------------------------------------------------------------
analyze_sleuth_subset_multi_D5('Day','==','D5','Treatment','D5_Treatment')
analyze_sleuth_subset_multi_D6('Day','==','D6','Treatment','D6_Treatment')
analyze_sleuth_subset_multi_D9('Day','==','D9','Treatment','D9_Treatment')
analyze_sleuth_subset_multi_D6_vs_D9('Day','==','D9','Treatment','D6_D9_Olap')
analyze_sleuth_subset_olaparib('Day','==','D9','Treatment','Olap_DMSO')
analyze_sleuth_subset_D5_D9('Day','==','D9','Treatment','Olap_DMSO_D5_D9')
