#!/bin/bash

#Author: Nathan T. Johnson
#Email: johnsonnathant@gmail.com, nathan_johnson@hms.havard.edu
#Purpose: Download pipeline environments to run on local machine
#Assumption: In main directory where 'local' directory is one level below
#To Run: ./install.sh

echo 'Setting Up CyCif Pipeline for Local Machine Run'

#setup directory (if not present make, otherwise ignore)
echo 'Make Folders'
mkdir -p environments
echo 'Download CyCif Pipeline Environment'
wget --output-document=cycif_pipeline.tar.gz https://www.dropbox.com/s/32sdgrcnoeji9ta/cycif_pipeline.tar.gz?dl=0
echo 'Download ImageJ Environment'
wget --output-document=ImageJ.tar.gz https://www.dropbox.com/s/gket7972dwhucss/ImageJ.tar.gz?dl=0
echo 'Downloading Ashlar Environment'
wget --output-document=ashlar.tar.gz https://www.dropbox.com/s/uhm7qrhvq5b6po6/ashlar.tar.gz?dl=0
echo 'Downloading cf25 Ashlar Environment'
wget --output-document=ashlar_cf25.tar.gz https://www.dropbox.com/s/iozzhqoqk6keqzj/ashlar_cf25.tar.gz?dl=0
echo 'Downloading Unet Environment'
wget --output-document=cycif-segment-tf-umap.tar.gz https://www.dropbox.com/s/qt4kmn0vvm6w007/unet.tar.gz?dl=0 
echo 'Downloading Segmenter via Clarence'
wget --output-document=segmenter.tar.gz https://www.dropbox.com/s/w9fniau7od8c6iv/segmenter.tar.gz?dl=0
echo 'Downloading Feature Extractor via HistoCat'
wget --output-document=histoCAT.tar.gz https://www.dropbox.com/s/8iawjlcwa9jo14o/histoCAT.tar.gz?dl=0
echo 'Uncompressing ImageJ'
tar -zxf ImageJ.tar.gz --directory ./environments
echo 'Uncompressing Ashlar'
tar -zxf ashlar.tar.gz --directory ./environments
echo 'Uncompressing cf25 Ashlar'
tar -zxf ashlar_cf25.tar.gz --directory ./environments
echo 'Uncompressing Unet'
tar -zxf cycif-segment-tf-umap.tar.gz --directory ./environments
echo 'Uncompressing Segmenter'
tar -zxf segmenter.tar.gz --directory ./environments
echo 'Uncompressing Feature Extractor'
tar -zxf histoCAT.tar.gz --directory ./environments
echo 'Uncompressing CyCif Pipeline'
tar -zxf ./cycif_pipeline.tar.gz --directory ./environments
echo 'Clean Up'
rm *.tar.gz
