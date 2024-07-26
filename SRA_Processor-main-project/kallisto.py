'''
Program to process Download, SRA extract, Delete SRA, compress Fastq files, run kallisto
to run: python sra_kallisto.py [file of SRR IDs to process]
#To add: ability to take in kallisto index, and path information

'''

import os
import sys
import pandas as pd

#path
path="/home/bionerd/Dropbox/Company/@Verne_Bioanalytics/@Assay_Development/"

#file
#file = pd.read_csv(''.join([path,'Data/RNA_Seq/SRA_Ids.csv']),names=['Sample'])
file = pd.read_csv(sys.argv[1],names=['Sample'])

for i in file.Sample:
    print("Processing:",i)
    #download data
    #cmd = ''.join(["prefetch -O ",path,"Data/RNA_Seq/SRA_Data/ ", i])
    #print(cmd)
    #os.system(cmd)

    #extract SRA data
    #cmd = ''.join([path,"Data/RNA_Seq/SRA_Data/",i])
    #print(cmd)
    #os.chdir(cmd)

    #cmd = ''.join(["fasterq-dump -S ",i])
    #print(cmd)
    #os.system(cmd)

    #delete SRA
    #cmd = "rm *.sra"
    #print(cmd)
    #os.system(cmd)

    #compress data
    #cmd = "gzip *.fastq"
    #print(cmd)
    #os.system(cmd)

    #Run Kallisto

    index =''.join([path,"Data/Cannabis_sativa_NCBI_sequences/Kallisto/Cannabis_V1.kallisto_index"])
    #gtf = ''.join([path,"/Reference/Homo_sapiens.GRCh38.101.gtf.gz"])
    output = ''.join([path,"Data/RNA_Seq/SRA_Kallisto_Data/",i])
    forward = ''.join([path,"Data/RNA_Seq/SRA_Data/",i,"/",i,"_1.fastq.gz"])
    reverse = ''.join([path,"Data/RNA_Seq/SRA_Data/",i,"/",i,"_2.fastq.gz"])

    cmd = ''.join(["kallisto quant -i ",index," -b 10 -t 4 -o ", output," ", forward," ", reverse])
    print(cmd)
    os.system(cmd)
