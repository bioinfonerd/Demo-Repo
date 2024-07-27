


# Overview of Services 	
Provider shall provide services listed in the table below to Customer. 

| Product Name        | 10X Genomics single cell RNA seq                                                                                               |
| :------------------ | :----------------------------------------------------------------------------------------------------------------------------- |
| Sample Type         | Frozen cells                                                                                                                   |
| Sample Number       | 6                                                                                                                              |
| Species Name        | Mouse                                                                                                                          |
| Sample origin       | USA                                                                                                                            |
| Sequencing Platform | DNBseq                                                                                                                         |
| Read Length         | 100PE                                                                                                                          |
| Data output         | 50,000 raw reads per single cell; 5,000-6,000 single cell is expected to be generated for each sample                          |
| Bioinformatics      | No bioinformatic analysis                                                                                                      |
| Turnaround Time     | 35-40 business days                                                                                                            |
| Additional Comments | We recommend 50,000 raw reads per single cell. If we isolate 5,000 single cell for each sample, 50K * 5,000=250 M reads/sample |
|                     |                                                                                                                                |


## Fresh Sample viability Prior to Frozen
| sample name | viable cell (million/ml) | viability | vol (ml) | Frozen media |
|:---------- |:----------------------- |:-------- |:------- |:-----------|
| CT26-5-1    | 1.9                      | 82.50%    | 1        | CS10         |
| CT26-5-2    | 1.9                      | 82.50%    | 1        | CS10         |
| CT26-5-3    | 1.9                      | 82.50%    | 1        | CS10         |
| CT26-5-4    | 4                        | 83.80%    | 1        | CS10         |
| CT26-5-5    | 4                        | 83.80%    | 1        | CS10         |
| CT26-5-6    | 4                        | 83.80%    | 1        | CS10         |



## Sample ID Table & QC Results from BGI (cDNA) 
| Sample # | BGI Sample ID | Customer  | Sample Type | Sample Volume (ul) | BGI Concentration (ng/ul) | BGI Total Quantity (ng) |
| :------- | :------------ | :-------- | :---------- | :----------------- | :------------------------ | :---------------------- |
|          |               | Sample ID |             |                    |                           |                         |
| 1        | 19159-19-01   | CT26-5-1  | cDNA        | 40                 | 12.10                     | 484.0                   |
| 2        | 19159-19-02   | CT26-5-2  | cDNA        | 40                 | 15.60                     | 624.0                   |
| 3        | 19159-19-03   | CT26-5-3  | cDNA        | 40                 | 15.80                     | 632.0                   |
| 4        | 19159-19-04   | CT26-5-4  | cDNA        | 40                 | 6.68                      | 267.2                   |
| 5        | 19159-19-05   | CT26-5-5  | N/A         | N/A                | N/A                       | N/A                     |
| 6        | 19159-19-06   | CT26-5-6  | cDNA        | 40                 | 5.14                      | 205.6                   |
|          |               |           |             |                    |                           |                         |



## Sample ID Table & QC Results (Post dead cell removal)
| Sample # | BGI-ID      | Customer Sample ID | Sample Type | Sample Volume (ul) | Total Cell Count  (per ul) | Cell Viability | Number of cells loaded |
| :------- | :---------- | :----------------- | :---------- | :----------------- | :------------------------- | :------------- | :--------------------- |
| 1        | 19159-19-01 | CT26-5-1           | cell        | 150                | 478                        | 0.780          | 12000                  |
| 2        | 19159-19-02 | CT26-5-2           | cell        | 150                | 454                        | 0.784          | 12000                  |
| 3        | 19159-19-03 | CT26-5-3           | cell        | 150                | 490                        | 0.707          | 12000                  |
| 4        | 19159-19-04 | CT26-5-4           | cell        | 150                | 174                        | 0.738          | 7517                   |
| 5        | 19159-19-05 | CT26-5-5           | cell        | 150                | 173                        | 0.550          | N/A                    |
| 6        | 19159-19-06 | CT26-5-6           | cell        | 150                | 233                        | 0.760          | 10066                  |



```
From: Yaya Wang <yaya_wang@h3biomedicine.com>				
Sent: Monday, May 3, 2021 12:44 PM				
To: Zonghui Peng <zonghui.peng@bgi.com>				
Cc: Yonghong Xiao <yonghong_xiao@h3biomedicine.com>; Kuan-Chun Huang <kuan-chun_huang@h3biomedicine.com>; Yaoyu Wang <yaoyu_wang@h3biomedicine.com>				
Subject: Re: timing for scRNA-seq study samples				
				
外部邮件/External Mail				
				
Hi Zonghui,				
Below is the sample information. Should the cell viability not pass QC after thawing, please do a dead cell removal before loading to the 10x.				
Thank you.				
sample name	viable cell (million/ml)	viability	vol (ml)	Frozen media
CT26-5-1	1.9	82.50%	1	CS10
CT26-5-2	1.9	82.50%	1	CS10
CT26-5-3	1.9	82.50%	1	CS10
CT26-5-4	4	83.80%	1	CS10
CT26-5-5	4	83.80%	1	CS10
CT26-5-6	4	83.80%	1	CS10

Yaya			
```

```
何韵秋(Yunqiu He) <heyunqiu1@bgi.com>
Jun 4, 2021, 6:05 AM
to Kuan-Chun, Yaya, Yonghong, me, 涂碧梦(Bimeng, Zonghui, 刘嘉颖(Jiaying

Dear Kuan，

 Hope you are well today. Data is ready now. In total, we have completed scRNA-seq study for your 5 samples for F21FTSUSAT0272-LIBuepR and $21,544 will be charged. Here is the path.


Path: s3://h3bioinf-share-bgi/F21FTSUSAT0272_LIBuepR_4Jun2021

Thanks for choosing BGI and we look forward to working with you again.
Please feel free to let me know if any questions you have.
Best，

Yunqiu
----
Nathan Johnson
Jun 4, 2021, 2:10 PM
to 刘嘉颖(Jiaying, Yonghong, 何韵秋(Yunqiu, Kuan-Chun, Yaya, me, 涂碧梦(Bimeng, Zonghui

Thank you Yungqiu! I can confirm we have it. 
```
