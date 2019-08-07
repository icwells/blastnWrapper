'''Configuration class for blastnWrapper'''

from glob import glob
import gzip
from multiprocessing import cpu_count
import os
import re
from unixpath import *

class Config():

	def __init__(self, infile, threads):
		self.threads = self.__setThreads__(threads)
		self.bin = ""
		self.genome = ""
		self.database = ""
		self.fastas = {}
		self.fastqs = {}
		self.outdir = ""
		self.resdir = ""
		self.results = {}
		self.__setConfig__(infile)
		self.__checkOutdir__()

	def __setThreads__(self, n):
		# Prevent too many threads from being called
		if int(n) > cpu_count():
			n = str(cpu_count())
		return n

	def __fastqToFasta__(self, outdir, infile):
		# Converts input fastq to fasta, unzips if needed
		print(("\tConverting {} to fasta...").format(os.path.split(infile)[1]))
		score = re.compile(r"[0-9]|>|/|\?")
		outfile = os.path.join(outdir, getFileName(infile) + ".fasta")
		if ".gz" in infile:
			f = gzip.open(infile, "rt")
		else:
			f = open(infile, "r")
		with open(outfile, "w") as out:
			for line in f:
				if len(line) <= 2:
					# Replace with blank line
					out.write("\n")
				elif line[0] == "@":
					out.write(line.replace("@", ">"))
				elif not score.match(line):
					out.write(line)
		f.close()
		return outfile

	def __getFileDict__(self, indir, ext, ext2=None):
		# Reads in dictionary of paired read files
		rm = []
		ret = {}
		indir = checkDir(indir, False)
		files = glob(indir + "*{}*".format(ext))
		if ext2 is not None:
			files.append(indir + "*{}*".format(ext))
		for i in files:
			# Get sample names and add to dict
			name = getSampleName(i)
			if name not in ret.keys():
				ret[name] = ["", ""]
			# Add forward/reverse files in order
			if "_R1_" in i:
				ret[name][0] = i
			else:
				ret[name][1] = i
		for k in ret.keys():
			# Make sure there are two files per sample
			if ret[k][0] == "" or ret[k][1] == "":
				rm.append(k)
		for i in rm:
			del ret[i]
		return ret

	def __checkOutdir__(self):
		# Makes sure proper output directory structure exists
		self.outdir = checkDir(self.outdir, True)
		# Get results directory and files
		self.resdir = checkDir(os.path.join(self.outdir, "blastResults"), True)
		self.results = self.__getFileDict__(self.resdir, "outfmt6")
		# Get fasta directory and files
		fastas = checkDir(os.path.join(self.outdir, "fastas"), True)
		self.fastas = self.__getFileDict__(fastas, "fasta", "fa")
		for k in self.fastqs.keys():
			if k not in self.fastas.keys():
				self.fastas[k] = [self.__fastqToFasta__(fastas, self.fastqs[k][0]), self.__fastqToFasta__(fastas, self.fastqs[k][1])]

	def __setConfig__(self, infile):
		# Reads input from config file
		if not os.path.isfile(infile):
			print("\n\t[Error] Cannot find configuation file. Exiting.\n")
			quit()
		print("\n\tReading configuration file...")
		with open(infile, "r") as f:
			for line in f:
				line = line.split("=")
				key = line[0].strip()
				val = line[1].strip()
				if key == "blast bin directory":
					self.bin = checkDir(val, False)
				elif key == "reference genome":
					self.genome = checkFile(val)
				elif key == "blast database":
					self.database = val
				elif key == "input fastq files":
					self.fastqs = self.__getFileDict__(val, "fastq", "fq")
				elif key == "output directory":
					self.outdir = checkDir(val, True)
