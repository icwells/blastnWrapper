'''Wraps calls to blastn'''

from config import Config
import os
from datetime import datetime
from argparse import ArgumentParser
from unixpath import *

def blastSeqs(config):
	# Calls blastx on all input and blastn on all hits with e < 10^-5
	print("\tRunning blastn against nucleotide database...")
	for k in config.fastas.keys():
		if k not in config.results.keys():
			print(("\tCalling blastn on {}...").format(k))
			config.results[k] = ["", ""]
			for idx, i in enumerate(config.fastas[k]):
				outfile = os.path.join(config.resdir, GetFileName(i) + ".outfmt6")
				cmd = ("blastn -query {} -db {} -num_threads {} -max_target_seqs 1 -outfmt 6 -out {}").format(i, config.db, config.threads, outfile)
				res = runProc(cmd)
				if res == True:
					# Add to results dict
					config.results[k][idx] = outfile
				else:
					print(("\t[Error]: Blastn failure on {} R{}.").format(k, idx))
	return config

def makeDB(config):
	# Makes protein and nucleotide blast databases
	print("\tConstructing BLAST nucleotide database...")
	cmd = ("makeblastdb -in {} -parse_seqids -dbtype nucl").format(config.genome)
	res = runProc(cmd)
	if not res:
		print("\tError: Failed to build BLAST nucleotide database.")

def exportBin(config):
	# Exports blast bin directory to linux path
	unipath.runProc(("export PATH=$PATH:{}").format(config.bin))

def getArguments():
	# Returns parsed arguemts
	parser = ArgumentParser("Wraps calls to blastn")
	parser.add_argument("--makedb", default = False, action = "store_true", help = "Makes database for blasting with provided fasta file.")
	parser.add_argument("--blastn", default = False, action = "store_true", help = "Calls blastn on provided input files. Runs summary by default.")
	parser.add_argument("--summary", default = False, action = "store_true", help = "Produces summary of directory of ublast results.")
	parser.add_argument("-c", default = "../config.txt", help = "Path to configuration file.")
	parser.add_argument("-t", default = "1", help = "Number of threads to run blast.")
	parser.add_argument("-e", default = "0.00001", help = "Maximum expected value.")
	parser.add_argument("-v", default = False, action = "store_true", help = "Print version info.")
	return parser.parse_args()

def version():
	print("\n\tublaster v0.1 is a package for managing loacl alignments with the ublast program.")
	print("\n\tCopyright 2019 by Shawn Rupp, Maley Lab, Biodesign Institute, Arizona State University.")
	print("\tThis program comes with ABSOLUTELY NO WARRANTY.")
	print("\n\tThis is free software, and you are welcome to redistribute it under certain conditions.\n")
	quit()

def main():
	start = datetime.now()
	args = getArguments()
	if args.v:
		version()
	# Load config file and append blast to path
	conf = Config(args.c, args.t, args.e)
	exportBin(config)
	if args.makedb:
		makeDB(config)
	elif args.blastn:
		config = blastSeqs(config)
		#summarize(config)
	elif args.summary:
		summarize(config)
	else:
		print("\n[Error] Please supply a valid command.\n")
		quit()
	print(("\tTotal runtime: {}\n").format(datetime.now() - starttime))

if __name__ == "__main__":
	main()
