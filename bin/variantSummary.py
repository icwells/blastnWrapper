'''Compares blast output to variants files and produces summary of overlap'''

from collections import OrderedDict
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

	def __setVariant__(self, row):
		# Reads row of variant file into dict by id and chromosome
		pid = row[self.vhead["Patient"]]
		c = row[self.vhead["Chr"]]
		start = row[self.vhead["Start"]]
		end = row[self.vhead["End"]]
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
					idx.append(i.pid)
					s.append(i.hits)
					v.append(i.row)
					# Get match data
					ids, matches = i.getMatches()
					ind.extend(ids)
					m.extend(matches)
		summary = pandas.DataFrame(s, index = idx, columns = ["AmpliseqHits"])
		hits = pandas.DataFrame(m, index = ind, columns = ["Chr", "VariantStart", "VariantEnd", "QueryID", "QueryStart", "QueryEnd"])
		var = pandas.DataFrame(v, index = idx, columns = list(self.vhead.keys())[1:])
		return summary, hits, var

	def __writeVariants__(self):
		# Writes comparison results to file
		print("\tWriting results to file...")
		summary, hits, var = self.__getOutput__()
		with pandas.ExcelWriter(self.outfile) as writer:
			summary.to_excel(writer, sheet_name = "Summary")
			hits.to_excel(writer, sheet_name = "Blast Hits")
			var.to_excel(writer, sheet_name = "Variants")

	def __compareVariants__(self, name):
		# Determines if blast results overlap with variants
		for c in self.results.keys():
			if c in self.variants[name].keys():
				for v in self.variants[name][c]:
					for i in self.results[name][c]:
						# Detmerine if any blast hits contain variant locus
						if self.results[name][c].inRange(v.start, v.end):
							self.variants[name][c].add(i.pid, i.start, i.end)

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

	def __setBlastResults__(self, pid, infile):
		# Reads in infile as a dictionary stored by chromosome (each file is one sample)
		first = True
		with open(infile, "r") as f:
			for line in f:
				if first == True:
					delim = getDelim(line)
					first = False
				row = line.strip().split(delim)
				pas = self.__evaluateRows__(row)
				if pas == True:
					c = row[self.bhead["subjectid"]]
					start = row[self.bhead["sstart"]]
					end = row[self.bhead["send"]]
					if c not in self.results.keys():
						self.results[c] = {}
					self.variants[c].append(Variant(pid, c, start, end, row))

	def __compareResults__(self):
		# Compares all blast result files to variants
		for k in self.blast.keys():
			for i in self.blast[k]:
				# Clear dict for next file
				self.results.clear()
				name = getSampleName(i)
				if name in self.variants.keys():
					print(("\tComparing results from {}...").format(name))
					self.__setBlastResults__(name, i)
					self.__compareVariants__(name)
