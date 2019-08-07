'''Stores variant data in class instance'''

class Variant():

	def __init__(self, pid, c, start, end, row):
		self.id = pid
		self.chromosome = c
		self.start = float(start)
		self.end = float(end)
		self.row = row[1:]
		self.hits = 0
		self.matches = []

	def getMatches(self):
		# Returns list of query/subject matches
		ids = []
		ret = []
		stem = [self.chromosome, self.start, self.end]
		for i in self.matches:
			idx.append(self.id)
			ret.append(stem.extend(i))
		return ids, ret

	def add(self, pid, start, end):
		# Indexes hit counter
		self.hits += 1
		self.matches.append([pid, start, end])

	def inRange(self, start, end):
		# Returns true if start and end are inside self.start/end
		if self.start <= start and self.end >= end:
			return True
		else:
			return False
