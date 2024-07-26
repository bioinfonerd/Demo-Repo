#!/usr/bin/env python3
# coding: utf8

# ==============================================================================
# Title: SRA_Processor
# Description: Program Launch
# Author: Nathan T. Johnson
# Contact: johnsonnathant@gmail.com
# Update Date: 2021.01.16
# Language release: python 3.8.5
# ==============================================================================

'''
Program to process Download, SRA extract, Delete SRA, compress Fastq files
to run: python sra_process.py [path to configuration file]
python sra_process.py configuration.cfg

configuration file example is X

Parallel processing of SRA file management


'''

# <editor-fold desc="Library">
import os
import sys
import pandas as pd
import os.path
import glob
import logging
import configparser
import multiprocessing as mp



# </editor-fold>

# <editor-fold desc="Function">

# process parameter configuration file
def process_configuration_file(path='./configuration.cfg'):
    parameters = configparser.ConfigParser()
    parameters.read(path)
    return (parameters)

# display contents of the parameter file
def config_printing(parameters):
    # log
    logging.info('Contents of Configuration File')

    # log the sections digested from config file
    #logging.info('Sections in configuration file', parameters.sections())
    logging.info(f"Sections in configuration file: {parameters.sections()}")

    if not parameters.sections():
        logging.info('Parameter section is null, check configuration file or path')
        # stop program
        sys.exit()
    else:
        # for each section display the section contents
        for s in parameters.sections():
            logging.info('section %s', str(s))
            for key in parameters[s]:
                value = parameters[s][key]
                logging.info('key: %s %s', str(key), str(value))

# read in csv file that organizes SRA data
def SRA_file_process(parameters):
    df = pd.read_csv(''.join([os.path.abspath(parameters['path']['path_to_file']), '/', parameters['data']['file']]),
                     sep=',')
    return (df)

# clean up leftover SRA tmp files
def cleanup(path):
    # Check if SRA tmp files exist
    if glob.glob(''.join([path,'*.sra.*'])).__len__() != 0:
        logging.info(f'SRA Tmp Files Identifed in {path}, removing')
        # clean up any leftover SRA tmp files
        cmd = ''.join(['rm ',path,'*.sra.*'])
        os.system(cmd)

# process just the directory information for the project
def project_directory(file, parameters, lib_strategy):
    # library type selection not detected
    if file[file['LibraryStrategy'] == lib_strategy].empty:
        logging.error('No Library Type match %s, check configuration file under data,file_to_process',
                      lib_strategy)
    else:
        file = file[file['LibraryStrategy'] == lib_strategy]
        # make directory for each SRA project
        for i in file["SRAStudy"].unique():
            # directory to make
            cmd = os.path.join(parameters['path']['sra_directory'], lib_strategy, i)
            logging.info('Creating Directory: %s', cmd)
            # make full path of directory if doesn't exist
            if not os.path.exists(cmd):
                os.makedirs(cmd)
            else:
                logging.info('Checking for leftover SRA tmp files')
                cleanup(cmd)

# make directory
def directory(file, parameters):
    # determine whether there is a directory at location, if not create
    # if so, check if new directories should be added
    logging.info('Making Necessary Directories')
    # build full directory
    if parameters['data']['file_to_process'] == 'ALL':
        logging.info('Making a directory for all filetypes')
        for LibStrategy in file['LibraryStrategy'].unique():
            project_directory(file, parameters, LibStrategy)
    # select files to process based on request
    else:
        logging.info('Making a directory just for %s', parameters['data']['file_to_process'])
        project_directory(file, parameters, parameters['data']['file_to_process'])

# identify the total samples to process, make a list of 2d list with file name and path
def file_processing(file, parameters):
    logging.info('Identifying Files to Consider')
    if parameters['data']['file_to_process'] == 'ALL':

        # assign sample name as key
        keys = file['Run'].tolist()

        # build list from path
        values = []

        for i in keys:
            # verifying 1 result
            if file[file['Run'] == i].__len__() != 1:
                logging.info('Not single file: %s', i)
            # add path
            values.append(''.join([parameters['path']['sra_directory'], '/',
                                   file[file['Run'] == i]['LibraryStrategy'].iloc[0], '/',
                                   file[file['Run'] == i]['SRAStudy'].iloc[0]]))

        # using list comprehension to create 2d array within 1d array
        files_to_consider = [[keys[i],values[i]] for i in range(len(keys))]

        # Printing resultant dictionary
        logging.info("Files to consider are: %s", str(files_to_consider))

    # or focused on one particular library
    else:
        # get files for just the selected library
        # assign sample name as key
        keys = file[file['LibraryStrategy'] == parameters['data']['file_to_process']]['Run'].tolist()

        # build list from path
        values = []
        for i in keys:
            # verifying 1 result
            if file[file['Run'] == i].__len__() != 1:
                logging.info('Not single file: %s', i)
            # add path
            values.append(''.join([parameters['path']['sra_directory'], '/',
                                   file[file['Run'] == i]['LibraryStrategy'].iloc[0], '/',
                                   file[file['Run'] == i]['SRAStudy'].iloc[0]]))

        # using list comprehension to create 2d array within 1d array
        files_to_consider = [[keys[i],values[i]] for i in range(len(keys))]

        # Printing resultant dictionary
        logging.info("Files to consider are: %s", str(files_to_consider))
    return (files_to_consider)

# always check if files are processed in case program is to be stopped and rerun in the middle
# [TODO] does not consider if file was in the middle of execution
def update(files_to_consider):
    # 1d list for each file to track process
    files_to_download = []
    files_to_process = []
    files_to_compress = []
    files_completed = []

    logging.info('Identifying what files have been downloaded, extracted, and compressed')
    logging.info('Identifying Files to Download')

    # check if sample is found
    # file needs to be downloaded if no sample name (could be .sra, .fastq. or .fastq.gz representing different phases)
    for key in files_to_consider:
        if ((glob.glob(''.join([key[1], '/', key[0], '.*'])).__len__() != 0) | (glob.glob(''.join([key[1], '/', key[0], '_*'])).__len__() != 0)):
            #logging.info('Identified %s files has been downloaded', key[0])
            next
        else:
            files_to_download.append(key)

    logging.info('Identifying Files to Extract')
    # check if sra file found, if so add to list
    for key in files_to_consider:
        if os.path.exists(''.join([key[1], '/', key[0], '.sra'])):
            # add to list
            files_to_process.append(key)

    logging.info('Identifying which Files to compress')
    # check if extracted file found, if so add to list
    for key in files_to_consider:
        if os.path.exists(''.join([key[1], '/', key[0], '.fastq'])):
            # add to list
            files_to_compress.append(key)

    logging.info('Identifying which Files are completed')
    # check if extracted file found, if so add to list
    for key in files_to_consider:
        #check for single end files
        if os.path.exists(''.join([key[1], '/', key[0], '.fastq.gz'])):
            # add to list
            files_completed.append(key)
        #check for paired end files
        else:
            if (os.path.exists(''.join([key[1], '/', key[0], '_1.fastq.gz'])) & os.path.exists(
                    ''.join([key[1], '/', key[0], '_2.fastq.gz']))):
                # add to list
                files_completed.append(key)

    return (files_to_download, files_to_process, files_to_compress,files_completed)

#statistics on files to process
def file_statistics(files_to_consider,files_to_download, files_to_process, files_to_compress,files_completed):
    logging.info('File Statistics')
    logging.info(f'To Consider: {str(len(files_to_consider))}')
    logging.info(f'To Download: {str(len(files_to_download))}')
    logging.info(f'To Extract: {str(len(files_to_process))}')
    logging.info(f'To Compress: {str(len(files_to_compress))}')
    logging.info(f'Completed: {str(len(files_completed))}')

################
#Main Functions#
################
#makes use of the global variable: file (df of sra file info)

# download SRA file
def download(sample):
    logging.info(f'Downloading: {sample[0]}')
    #download
    cmd = ''.join(["prefetch -O ", sample[1]," --max-size ",
                   str(parameters['download']['max_size'])," ", sample[0]])
    print(cmd)
    os.system(cmd)

    #[TODO] need to implement whether download was successful, for now just print complete
    print(f'Downloading Completed: {sample[0]}')


# extract SRA file
def extract(sample):
    logging.info(f'Extracting: {sample[0]}')
    # change directory
    cmd = ''.join([sample[1]])
    print(cmd)
    os.chdir(cmd)

    #check if sample with SRA present
    if os.path.exists(''.join([sample[1], '/', sample[0], '.sra'])):
        # extract SRA data (newest version does not need to consider if output file should be single or paired end)
        cmd = ''.join(["fasterq-dump -e ",
                       str(parameters['fasterqdump']['threads']),
                       " -m ",str(parameters['fasterqdump']['memory']),
                       " --temp ", str(''.join([sample[1], '/', sample[0]])),
                       " -S ", sample[0],'.sra'])
        print(cmd)
        os.system(cmd)

        #clean up tmp directory
        cmd = ''.join(['rm -r ', str(''.join([sample[1], '/', sample[0]]))])
        print(cmd)
        os.system(cmd)

    else:
        logging.info(f'SRA File not present: {sample[0]}')

#delete SRA file & compress
def compress(sample):
    logging.info(f'Compressing: {sample[0]}')
    # change directory
    cmd = ''.join([sample[1]])
    print(cmd)
    os.chdir(cmd)

    #check if SRA file present & extracted fastq data
    #checking SRA file present
    if os.path.exists(''.join([sample[1], '/', sample[0], '.sra'])):
        # checking fastq data (check sample file whether its paired or single end files)
        if file[file['Run']==sample[0]]['LibraryLayout'].iloc[0] == 'PAIRED':
            # assuming only options are PAIRED and SINGLE
            # checking if both files extracted
            if (os.path.exists(''.join([sample[1], '/', sample[0], '_1.fastq'])) & os.path.exists(''.join([sample[1], '/', sample[0], '_2.fastq']))):

                # delete SRA
                cmd = ''.join(["rm ", sample[0], '.sra'])
                print(cmd)
                os.system(cmd)

                # compress data
                cmd = ''.join(["pigz -p ", str(parameters['compress']['threads'])," ",
                               sample[0], '*.fastq'])
                print(cmd)
                os.system(cmd)

        else:
            #assuming only options are PAIRED and SINGLE
            if os.path.exists(''.join([sample[1], '/', sample[0], '.fastq'])):
                # delete SRA
                cmd = ''.join(["rm ", sample[0], '.sra'])
                print(cmd)
                os.system(cmd)

                # compress data
                cmd = ''.join(["pigz -p ", str(parameters['compress']['threads'])," ",
                               sample[0], '*.fastq'])
                print(cmd)
                os.system(cmd)
            else:
                logging.info(f'Expected fastq file not present (assuming single read type): {sample[0]}')
    else
        logging.info(f'SRA file not present, perhaps stopped mid compressing?: {sample[0]}')
        logging.info(f'Checking if fastq file needs compressing: {sample[0]}')

        # checking fastq data (check sample file whether its paired or single end files)
        if file[file['Run']==sample[0]]['LibraryLayout'].iloc[0] == 'PAIRED':
            # assuming only options are PAIRED and SINGLE
            # checking if both files extracted
            if (os.path.exists(''.join([sample[1], '/', sample[0], '_1.fastq'])) & os.path.exists(
                    ''.join([sample[1], '/', sample[0], '_2.fastq']))):
                # compress data
                cmd = ''.join(["pigz -p ", str(parameters['compress']['threads']), " ",
                               sample[0], '*.fastq'])
                print(cmd)
                os.system(cmd)

        else:
            # assuming only options are PAIRED and SINGLE
            if os.path.exists(''.join([sample[1], '/', sample[0], '.fastq'])):
                # compress data
                cmd = ''.join(["pigz -p ", str(parameters['compress']['threads']), " ",
                               sample[0], '*.fastq'])
                print(cmd)
                os.system(cmd)
            else:
                logging.info(f'Expected fastq file not present (assuming single read type): {sample[0]}')

#####################################
#Combine Functions based on Use Case#
#####################################
#makes use of the global variable: file (df of sra file info)

#extract, and compress
def extract_compress(sample):
    extract(sample)

    compress(sample)
# </editor-fold>

if __name__ == '__main__':
    # change directory
    # os.chdir('/home/bionerd/Dropbox/Company/@Verne_Bioanalytics/git/SRA_Processor')
    # logging output file
    # overwrite log each run of program
    # TODO log file will have an appended time so new log each run, currently just overwrites previous log run
    logging.basicConfig(  # ALL, DEBUG, INFO, ERROR, FATAL
        filename='log.txt',
        filemode='w',
        level=logging.INFO)

    # process parameter file
    parameters=process_configuration_file(path=sys.argv[1])
    # parameters = process_configuration_file(path='./configuration.cfg')

    # log parameters files
    logging.info('Configuration file processed')

    # what parameters digested
    config_printing(parameters)

    # read in SRA file
    file = SRA_file_process(parameters)

    # prep directory structure
    directory(file, parameters)

    # what files should be considered
    files_to_consider = file_processing(file, parameters)

    # maintain 3 separate pools to download, extract, and compress
    # Due to its not wise to communicate between pools, the process is:
    # 1) all files to download go through download, extract, and compress
    # 2) then repeat for extract and compress, allows balance between able to repeat in middle of program, but maintain concurrency

    ##########
    #Download#
    ##########

    # identify files to download, extract, and compress
    files_to_download, files_to_process, files_to_compress, files_completed = update(files_to_consider)

    #statistics to process
    file_statistics(files_to_consider, files_to_download, files_to_process, files_to_compress, files_completed)
    
    #skip process if already downloaded

    if len(files_to_download) == 0:
        logging.info('All Files Downloaded')
    else:
        # download thread pool
        logging.info('Starting Download Process')
        with mp.Pool(processes=int(parameters['download']['download_parallel_files'])) as pool:  # Parallelize
            pool.map(download, files_to_download)

    #########
    #Extract#
    #########
    logging.info('Checking For Files to Extract and Compress Files')
    # identify files to download, extract, and compress
    files_to_download, files_to_process, files_to_compress, files_completed = update(files_to_consider)

    #statistics to process
    file_statistics(files_to_consider, files_to_download, files_to_process, files_to_compress, files_completed)

    # Extract & Compress Process pool
    logging.info('Starting Extract & Compress Process')

    if len(files_to_process) == 0:
        logging.info('All Files Extracted')
    else:
        # extract, compress thread pool
        logging.info('Starting Extract and Compress Process')
        with mp.Pool(processes=int(parameters['fasterqdump']['fasterqdump_parallel_files'])) as pool:  # Parallelize
            pool.map(extract_compress, files_to_process)

    ###########
    #Compresss#
    ###########
    logging.info('Checking Any Leftover Compressed Files')
    #in case any other files needing to be compressed

    # identify files to download, extract, and compress
    files_to_download, files_to_process, files_to_compress, files_completed = update(files_to_consider)

    #statistics to process
    file_statistics(files_to_consider, files_to_download, files_to_process, files_to_compress, files_completed)

    if len(files_to_compress):
        logging.info('All Files Extracted')
    else:
        # Compress any leftover files Process pool
        for i in files_to_compress:
            compress(i)
    
    #######
    #Final#
    #######
    logging.info('Final Statistics')
    # identify files to download, extract, and compress
    files_to_download, files_to_process, files_to_compress, files_completed = update(files_to_consider)
    #statistics to process
    file_statistics(files_to_consider, files_to_download, files_to_process, files_to_compress, files_completed)
