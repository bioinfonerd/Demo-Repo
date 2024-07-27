# [company name]-scRNAseq-Analysis
### Finalized 2023.09.01 by Nathan Johnson

### Analyze scRNA-Seq data for [company name] with the following experiment key. 

- Client wanted data analyzed as a whole cohort plus additional information from "working-doc-senescence.pdf"
- Analysis Files Assumed Data processed via Cell Ranger
- Originally analysis was done with a master R script, but with the additional requests of additional datasets, the process was split into 4 R scripts that would iterate through the different portions (Seurat Object formation, Normalization & Highest Variable Genes, UMAP, Clustering, Biomarker identification, and Cell Type Calling and Annotation

#### Experiment Key
Condition Key Senescence Study 
Experiment SE2: 
Harvest day at 28 post bleomycin (on March 17, 2022) 
Material for sequencing: Right lobe from 1 mouse/group. Shipping frozen in cryomedium on dry ice in minced form 

•              1X Control + vehicle = “SE2 control” (SE2_control_S2_L002) [CONTROL] 

•              1X Bleo + vehicle = “SE2 bleo” (SE2_bleomycin_S3_L002) [DISEASE CONTROL] 

•              1X Bleo + BCLXL inhibitor = “SE2 Drug 1” (SE2_Drug1_S4_L003) [DISEASE + TREATMENT] 

•              1X Bleo + MCL1 inhibitor  = “SE2 Drug 2” (SE2_Drug2_S5_L003) [DISEASE + MCL1 TREATMENT] 

•              1X Bleo + EphrinB2 antibody = “SE2 Drug 3A” (SE2_Drug_3A_S6_L003) [DISEASE +  EPHRIN TREATMENT] 
 

Experiment SE3: 
Harvest day at 16 post bleomycin (on April 21, 2022) 
Material for sequencing: Right lobes from 2-3 mice/group (pooled samples).  Shipping frozen in cryomedium on dry ice in minced form 
•              2X Controls  = “SE3 control” (SE3_control_S7_L003) [CONTROL] 

•              2X Bleo = “SE3 bleo” (SE3_bleo_S8_L003) [DISEASE CONTROL] 

•              3X Bleo + BCL-XL- inhibitor = “SE3 Drug A” (SE3_Drug_A_S9_L003) [DISEASE + BCL-XL TREATMENT] 

•              3X Bleo + MCL-1 inhibitor = “SE3 Drug B” (SE3_Drug_B_S10_L004) [DISEASE + MCL-1 TREATMENT] 

