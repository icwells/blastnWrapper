'''Wraps calls to blastn'''

from argparse import ArgumentParser
from config import Config
from datetime import datetime
import os
import re
from unixpath import *
from variantSummary import VariantSummary

def getSampleName(infile):
	# Returns sample name from filename
	name = getFilename(infile)
	name = name.split("_")[0]
	name = name.split("-")
	# Add DCIS and sample number
	ret = name[0] + name[1]
	if re.match(r"\d", name[2]) is not None:
		# Add alpha-numeric suffix
		ret = ("{}_{}").format(ret, name[2])
	return ret

def blastSeqs(config):
	# Calls blastx on all input and blastn on all hits with e < 10^-5
	print("\tRunning blastn against nucleotide database...")
	for k in config.fastas.keys():
		if k not in config.results.keys():
			print(("\tCalling blastn on {}...").format(k))
			config.results[k] = ["", ""]
			for idx, i in enumerate(config.fastas[k]):
				outfile = os.path.join(config.resdir, getFileName(i) + ".outfmt6")
				cmd = ("{}blastn -query {} -db {} -num_threads {} -max_target_seqs 1 -outfmt 6 -out {}").format(config.bin, i, config.database, config.threads, outfile)
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
	cmd = ("{}makeblastdb -in {} -parse_seqids -dbtype nucl").format(config.bin, config.genome)
	res = runProc(cmd)
	if not res:
		print("\tError: Failed to build BLAST nucleotide database.")

def getArguments():
	# Returns parsed arguemts
	parser = ArgumentParser("Wraps calls to blastn")
	parser.add_argument("--makedb", default = False, action = "store_true", help = "Makes database for blasting with provided fasta file.")
	parser.add_argument("--blastn", default = False, action = "store_true", help = "Calls blastn on provided input files. Runs summary by default.")
	parser.add_argument("--summary", default = False, action = "store_true", help = "Produces summary of directory of ublast results.")
	parser.add_argument("-c", default = "../config.txt", help = "Path to configuration file.")
	parser.add_argument("-t", default = "1", help = "Number of threads to run blast.")
	parser.add_argument("-i", help = "Path to input variants file.")
	parser.add_argument("-p", type = float, default = 95.0, help = "Minimum percent identity.")
	parser.add_argument("-e", type = float, default = 0.00001, help = "Maximum expected value.")
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
	config = Config(args.c, args.t)
	if args.makedb:
		makeDB(config)
	elif args.blastn:
		config = blastSeqs(config)
		if args.i:
			VariantSummary(args.i, config.outdir, args.p, args.e, config.results)
	elif args.summary and args.i:
		VariantSummary(args.i, config.outdir, args.p, args.e, config.results)
	else:
		print("\n[Error] Please supply a valid command.\n")
		quit()
	print(("\tTotal runtime: {}\n").format(datetime.now() - start))

if __name__ == "__main__":
	main()
