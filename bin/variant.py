'''Stores variant data in class instance'''

class Variant():

	def __init__(self, pid, c, start, end, row, name=None):
		self.id = pid
		self.chromosome = c
		self.name = name
		self.start = float(start)
		self.end = float(end)
		self.row = row[1:]
		self.matches = []

	def getMatches(self):
		# Returns list of query/subject matches
		ids = []
		ret = []
		for i in self.matches:
			stem = [self.chromosome, self.start, self.end]
			stem.extend(i)
			ids.append(self.id)
			ret.append(stem)
		return ids, ret

	def add(self, pid, start, end):
		# Appends matches
		self.matches.append([pid, start, end])

	def inRange(self, start, end):
		# Returns true if start and end are inside self.start/end
		if self.start <= start and self.end >= end:
			return True
		else:
			return False
