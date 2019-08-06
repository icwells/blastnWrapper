'''Stores variant data in class instance'''

class Variant():

	def __init__(self, row):
		self.row = row
		self.hits = 0
		self.matches = []

	def __str__(self):
		# Returns formatted string
		return ("{},{},{}").format(",".join(self.row), self.hits, ",".join(self.matches))

	def add(self, hit):
		# Appends hit to list and indexes hit counter
		self.matches.append(hit)
		self.hits += 1
		
