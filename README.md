# blastnWrapper is a Python 3 wrapper for the blastn function of ncbi blast+.  
These scripts are still under developement, but are available for use.  

Copyright 2019 by Shawn Rupp

1. [Dependencies](#Dependencies)  
2. [Installation](#Installation)  
3. [Usage](#Usage)  

## Dependencies:  

### Blast+  
Download NCBI Blast+ [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)  

## Installation  

	git clone https://github.com/icwells/blastnWrapper.git  

## Usage  

### Config File  
Make a copy of example_config.txt and fill in the fields with the path to appropriate file or directory.  

	blast bin directory		Directory containing Blast binaries.  
	reference genome		Path to Reference genome used for construcitng database.  
	blast database			Path and name of blast database (do not include a file extension).  
	input fastq files		Directory containing query of fastq files. They may or may not be gzipped.  
	output directory		Path to directory where output will be written.  

### Commands  

	python blastnWrapper.py --{makedb/blastn/summary} -c path_to_config

	-h, --help		show this help message and exit.  
	--makedb		Makes database for blasting with provided fasta file. Output is written to the blast database given in the config file.  
	--blastn		Calls blastn on provided input files. Runs summary by default.  
	--summary		Produces summary of directory of ublast results.  
	-c C			Path to configuration file.  
	-t 1			Number of threads to run blast.  
	-i I			Path to input variants file.  
	-p 0.95			Minimum percent identity.  
	-e 0.00001		Maximum expected value.  
	-v				Print version info.  

