# BCMA-SMLA DFCI CRISPR Data Analysis

## Overview of the data analysis
- CRISPR-ko using the [Brunello library](https://www.addgene.org/pooled-library/broadgpp-human-knockout-brunello/) covering the whole transcriptome.
- Used [KMS11](https://depmap.org/portal/cell_line/ACH-000714) cell line
- Resulting libraries were sequenced on NextSeq500
- Fastq files analyzed using MAGeCK workflow to obtain read counts and using MAGeCK MLE to determine essential genes.
- Contrasts performed are:
  
|control|condition|
|-|-|
HJFK7BGXK_KMS11_Brunello_9w_belamaf | HJFK7BGXK_KMS11_Brunello_9w_IgG1
HJFK7BGXK_KMS11_Brunello_9w_payload	| HJFK7BGXK_KMS11_Brunello_9w_DMSO
HJFK7BGXK_KMS11_Brunello_9w_SMLA7503 |	HJFK7BGXK_KMS11_Brunello_9w_IgG1
HJFK7BGXK_KMS11_Brunello_9w_SMLA7503 |	HJFK7BGXK_KMS11_Brunello_9w_belamaf
HJFK7BGXK_KMS11_Brunello_9w_belamaf	| HLFYLBGXY_BRUN_PDNA
HJFK7BGXK_KMS11_Brunello_9w_payload	| HLFYLBGXY_BRUN_PDNA
HJFK7BGXK_KMS11_Brunello_9w_SMLA7503 |	HLFYLBGXY_BRUN_PDNA

- Data files are availble in box folder [H3-Mitsiades SRA (H3B-7503)](https://h3biomedicine.box.com/s/epc2n9pqtglfkujdjrikvnz75s5xbtui)/
- GSEA was performed for all contrasts using **1)** C2; **2)** GO; **3)** all MSigDB


## Key questions
1. How is H3B-7503 different from Belamaf?
2. What are the H3B-7503 killing and resistance mechanism
3. What are the difference between SMLA and payload?

## What are the main H3B-7503 resistance pathways
### SMLA vs IgG
#### Positive Selection (H3B-7503 resistance)
Enriched canonical pathways showing gene overlaps between sets:
![H3B-7503 resistance CP groups](Figures/CRISPR_H3B7503vsIgG_pos_top200_overlap_path.png)

Key genes/complexes:
- **Group 1:**
  - NF$\kappa$B IA/IB
  - PTEN
  - TNFR1/2
  - GSK3$\beta$
  - TNFRSF17
  - [BIRC3](https://cellandbioscience.biomedcentral.com/articles/10.1186/s13578-020-00521-0)
    - "BIRC3 genetic inactivation due to deletions or point mutations is consistently associated to shorter progression free survival and poor prognosis in chronic lymphocytic leukemia patients. BIRC3 inactivation has also been associated to chemoimmunotherapy resistance."
  - cIAP
- **Group 2:**
  - EIF3 ("EIF3B", "EIF3C", "EIF3F", "EIF3H", "EIF3M")
  - EIF4 ("EIF4A1","EIF4A2")
  - RPS2
  - GSK3B
- **Group 3:**
  - Ferroptosis Singaling Pathway:
    - CASL4, BECN1, G3BP1, H2AC21, H2AZ,SQSTM1
  - Stearate Biosynthesis 1:
    - ACSL4
    - LIPT2
    - SRD5A3
  - 
- **Group 4:**
    - GPI
    - MPI



- [**Dolichyl-diphosphooligosaccharide Biosynthesis**](https://reports.ingenuity.com/rs/report/cpathway?id=ING%3A7l6h7)
  - ALG3
  - DPM1
- [**Putrescine Biosynthesis III**](https://reports.ingenuity.com/rs/report/cpathway?id=ING%3A7l6mj)
  -ODC1

#### Negative Selection (H3B-7503 sensitive)
Enriched canonical pathways showing gene overlaps between sets:
![H3B-7503 sensitive CP groups](Figures/CRISPR_H3B7503vsIgG_neg_top200_overlap_path.png)

Key genes/complexes:
- **Group 1:**
  - $\gamma$-secretase(APH1A, PSEN1, PSEN2, PSENEN, NCSTN)
  - CD3D
  - CAPN10 (Calpain)
  - CHP1
  - BIRC5

- **Group 2:**  
  - MLCP (myosin light chain phosphotase)(PPP1CB, PPP1R12A)
  - LIMK2
  - WAVE (WASF1)
  - NMDAR(GRIA3A, GRIA4)

- **DNA Methylation and Transcriptional Repression Signaling:**
  - DNMT3A
  - RBBP7

- **Glutathione Biosynthesis:**
  - GCLM

- **2-Ketoglutarate Dehydrogenase Complex:**
  - DLST




### Belamaf vs IgG
#### Positive Selection (Belamaf resistance)
![Belamaf vs IgG pos](Figures/CRISPR_BelamafvsIgG_pos_top200_overlap_path.png)


#### Negative Selection (Belamaf sensitive)


## What are the resistance mechanism
- BCL2 pathway
- 


## How is H3B-7503 different from Belamaf?


- KMS11 has high BCMA level compare to OPM2
- KMS11 sensive but reduced sensitive vs OPM2



## Expand the number of cell line models

- Check PoC slide 5 ![slide 5](Figures/20211209_PoC_slide5.png)
- Check PoC slide 7 ![slide 7](Figures/20211209_PoC_slide7.png)
  - Cell lines used here are less responsive to H3B-7503 and Belamaf in experiments done in Crown versus in-house.
  - MM1.R cell line was not responsive
  - Hypothesis: overall faster tumor growth in studies at Crown may blunt the actvitiy
- Check PoC slide 8 ![slide 8](Figures/20211209_PoC_slide8.png)
- Check PoC slide 10 ![slide 10](Figures/20211209_PoC_slide10.png)
  - Additional resistant cell lines (SKMM2, KMM1, HuNS1, KMS34)
  - 

