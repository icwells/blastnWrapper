'''Compares blast output to variants files and produces summary of overlap'''

from collections import OrderedDict
from config import getSampleName
import os
import pandas
from unixpath import *
from variant import Variant

class VariantSummary():

	def __init__(self, infile, outdir, pid, evalue, blastfiles):
		self.infile = infile
		self.outfile = os.path.join(outdir, "variantsBlastSummary.xlsx")
		self.pid = pid
		self.e = evalue
		self.vhead = OrderedDict()
		self.variants = {}
		self.bhead = {"queryid":0, "subjectid":1, "pidentity":2, "length":3, "mismatches":4, "gapopenings":5, 
					"qstart":6, "qend":7, "sstart":8, "send":9, "evalue":10, "bitscore":11}
		self.blast = blastfiles
		self.results = {}
		self.__setVariants__()
		self.__compareResults__()
		self.__writeVariants__()

	def __setHeader__(self, row):
		# Returns header dict for variant files
		for idx, i in enumerate(row):
			self.vhead[i.strip()] = idx

	def __setChromosome__(self, val):
		# Removes decimal from chromosome number
		if ".0" in val:
			val = val.split(".")[0]
		return val

	def __setVariant__(self, row):
		# Reads row of variant file into dict by id and chromosome
		pid = row[self.vhead["Patient"]]
		c = self.__setChromosome__(row[self.vhead["Chr"]])
		start = row[self.vhead["Start"]]
		end = row[self.vhead["End"]]
		name = row[self.vhead["Name"]]
		if pid not in self.variants.keys():
			self.variants[pid] = {}
		if c not in self.variants[pid].keys():
			self.variants[pid][c] = []
		self.variants[pid][c].append(Variant(pid, c, start, end, row))

	def __setVariants__(self):
		# Reads in dict of variants from infile
		checkFile(self.infile)
		print("\n\tReading variants file...")
		first = True
		with open(self.infile, "r") as f:
			for line in f:
				line = line.strip()
				if first == False:
					self.__setVariant__(line.split(delim))
				else:
					delim = getDelim(line)
					self.__setHeader__(line.split(delim))
					first = False

#-----------------------------------------------------------------------------

	def __getOutput__(self):
		# Returns data frame of summary data, alignment data, and input variants file
		s = []
		m = []
		v = []
		idx = []
		ind = []
		for k in self.variants.keys():
			for c in self.variants[k].keys():
				for i in self.variants[k][c]:
					# Store ids as seperate index for pandas
					idx.append(i.id)
					s.append([c, i.start, i.end, i.name, len(i.matches)])
					v.append(i.row)
					# Get match data
					ids, matches = i.getMatches()
					ind.extend(ids)
					m.extend(matches)
		summary = pandas.DataFrame(s, index = idx, columns = ["Chr", "Start", "End", "Name", "Coverage"])
		hits = pandas.DataFrame(m, index = ind, columns = ["Chr", "VariantStart", "VariantEnd", "QueryID", "QueryStart", "QueryEnd"])
		var = pandas.DataFrame(v, index = idx, columns = list(self.vhead.keys())[1:])
		return summary, hits, var

	def __writeVariants__(self):
		# Writes comparison results to file
		print("\n\tWriting results to file...")
		summary, hits, var = self.__getOutput__()
		with pandas.ExcelWriter(self.outfile) as writer:
			summary.to_excel(writer, sheet_name = "Summary")
			hits.to_excel(writer, sheet_name = "Blast Hits")
			var.to_excel(writer, sheet_name = "Variants")

	def __compareVariants__(self, name):
		# Determines if blast results overlap with variants for given sample
		for c in self.results.keys():
			if c in self.variants[name].keys():
				for idx, v in enumerate(self.variants[name][c]):
					for i in self.results[c]:
						# Detmerine if any blast hits contain variant locus
						if i.inRange(v.start, v.end):
							self.variants[name][c][idx].add(i.id, i.start, i.end)

	def __evaluateRows__(self, row):
		# Returns True if pid and evalue pass
		ret = False
		try:
			pid = float(row[self.bhead["pidentity"]])
			e = float(row[self.bhead["evalue"]])
			if pid >= self.pid and e < self.e:
				ret = True
		except TypeError:
			pass
		return ret

	def __setBlastResults__(self, name, infile):
		# Reads in infile as a dictionary stored by chromosome (each file is one sample)
		first = True
		with open(infile, "r") as f:
			for line in f:
				if first == True:
					delim = getDelim(line)
					first = False
				row = line.strip().split(delim)
				c = row[self.bhead["subjectid"]]
				pas = self.__evaluateRows__(row)
				if pas == True and c in self.variants[name].keys():
					# Only proceed if there is sufficient match quality and chromosome is present in variants
					qid = row[self.bhead["queryid"]]
					start = row[self.bhead["sstart"]]
					end = row[self.bhead["send"]]
					if c not in self.results.keys():
						self.results[c] = []
					self.results[c].append(Variant(qid, c, start, end, row))

	def __getSampleID__(self, filename):
		# Attempts to resolve result of getSampleName with variants keys
		ret = getSampleName(filename)
		if ret in self.variants.keys():
			return ret
		elif "_" in ret:
			ret = ret.split("_")[0]
			if ret in self.variants.keys():
				return ret
		return None			

	def __compareResults__(self):
		# Compares all blast result files to variants
		for k in self.blast.keys():
			name = self.__getSampleID__(self.blast[k][0])
			if name is not None:
				for idx, i in enumerate(self.blast[k]):
					# Clear dict for next file
					self.results.clear()
					print(("\tComparing results from {} R{}...").format(name, idx+1))
					self.__setBlastResults__(name, i)
					self.__compareVariants__(name)
