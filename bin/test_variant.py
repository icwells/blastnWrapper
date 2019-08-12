'''Tests variant and variantSummary classes'''

import pytest
import unixpath
from variant import Variant

def getVariants():
	# Returns initialized variants for testing
	ret = {}
	pid = "DCIS_1"
	rows = ["", ""]
	ret["1"] = [Variant(pid, "1", "100.0", "200.0", rows)]
	ret["1"].append(Variant(pid, "1", "1025", "1119", rows))
	ret["2"] = [Variant(pid, "2", "25006", "25124", rows)]
	ret["X"] = [Variant(pid, "X", "90045", "90157.5", rows)]
	return ret

def test_inRange():
	# Tests variant.inRange
	variants = getVariants()
	cases = {"1": [155, 156, True], "2": [8875, 8875, False], "X": [90065, 90065, True], "Y": [1500, 1506, False]}
	for k in cases.keys():
		if k in variants.keys():
			actual = False
			for i in variants[k]:
				res = i.inRange(cases[k][0], cases[k][1])
				if res == True:
					actual = res
			assert actual == cases[k][2]
		else:
			assert cases[k][2] == False
