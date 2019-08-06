'''Compares blast output to variants files and produces summary of overlap'''

from collections import OrderedDict
import os
from unixpath import *
from variant import Variant

class VariantSummary():

	def __init__(self, infile, outdir, pid, evalue, blastfiles):
		self.infile = infile
		self.outfile = os.path.join(outdir, "variantsBlastSummary.csv")
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
		# Reads row of variant file into dict
		c = row[self.vhead["Chr"]]
		start = float(row[self.vhead["Start"]])
		end = float(row[self.vhead["End"]])
		if c not in self.variants.keys():
			self.variants[c] = {}
		if start not in self.variants[c].keys():
			self.variants[c][start] = {}
		self.variants[c][start][end] = Variant(row)

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

	def __writeVariants__(self):
		# Writes comparison results to file
		total = 0
		hits = 0
		head = ("{},#Hits,HitIDs\n").format(",".join(self.vhead.keys()))
		print("\tWriting results to file...")
		with open(self.outfile, "w") as out:
			out.write(head)
			for c in self.variants.keys():
				for start in self.variants[c].keys():
					for end in self.variants[c][start].keys():
						out.write(str(self.variants[c][start][end]) + "\n")
						total += 1
						if self.variants[c][start][end].hits > 0:
							hits += 1
		print(("\tFound {} matches for {} total variants.").format(hits, total))

	def __binarySearch__(self, loci, target, less):
		# Performs binary search and returns closest lesser match for start and greater match for end
		strt = 0
		end = len(loci) - 1
		while strt <= end:
			md = int(start+end/2)
			val = loci[md]
			if target == val:
				return val
			elif less == True and target > val:
				end = md + 1
			elif less == True and target < val:
				start = md - 1
			elif target > val:
				start = md + 1
			else:
				end = md -1
		# Return closest match
		return val

	def __getLocus__(self, c, s, e):
		# Returns best overlapping locus
		start = None
		end = None
		if s in self.variants[c].keys():
			start = s
		else:
			start = self.__binarySearch__(list(self.variants[c].keys()), s, True)
		if e in self.variants[c][start].keys():
			end = e
		else:		
			end = self.__binarySearch__(list(self.variants[c][start].keys()), e, False)
		return start, end

	def __compareVariants__(self, name):
		# Determines if blast results overlap with variants
		for c in self.results.keys():
			if c in self.variants.keys():
				for s in self.results[c].keys():
					for e in self.results[c][s].keys():
						print(c, s, e)
						# Get best match and add to variant
						start, end = self.__getLocus__(c, s, e)
						if start is not None and end is not None:
							qid = self.results[c][s][e][self.bhead["queryid"]]
							self.results[c][start][end].add(("{}:{}").format(name, qid))

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
		# Reads in infile as a dictionary
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
					start = float(row[self.bhead["sstart"]])
					end = float(row[self.bhead["send"]])
					if c not in self.results.keys():
						self.results[c] = {}
					if start not in self.results[c].keys():
						self.results[c][start] = {}
					self.results[c][start][end] = row	

	def __compareResults__(self):
		# Compares all blast result files to variants
		for k in self.blast.keys():
			for i in self.blast[k]:
				# Clear dict for next file
				self.results.clear()
				name = getFileName(i)
				print(("\tComparing results from {}...").format(name))
				self.__setBlastResults__(name, i)
				self.__compareVariants__(name)
