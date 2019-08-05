'''This script contains common functions for manipulating paths in unix platforms.'''

import os
from sys import stderr
from subprocess import Popen
from shlex import split
from sys import stderr, stdout

def __printError__(msg):
	# Prints formatted error message and exits
	print(("\n\t[Error] {}. Exiting.\n").format(msg), file = stderr)
	quit()

def checkDir(path, make = False):
	# Checks if path exists, throws error if make is false, makes dir if true, returns formatted name
	if type(path) != str:
		__printError__(("{} is not a string").format(infile))
	if path[-1] != os.path.sep:
		path += os.path.sep
	if os.path.isdir(path) == False:
		if make == False:
			__printError__(("{} not found").format(path))
		else:
			os.mkdir(path)
	return path

def checkFile(infile):
	# Raises if infile is not found
	if type(infile) != str:
		__printError__(("{} is not a string").format(infile))
	if not os.path.isfile(infile):
		__printError__(("{} not found").format(infile))

def getExt(name):
	# Returns extension from filename
	return name[name.rfind(".")+1:]

def getFileName(path):
	# Returns filename sans path and extension
	return path[path.rfind(os.path.sep)+1:path.find(".")]

def getParent(path):
	# Returns parent directory from path
	if "." in path and path[-1] != os.path.sep:
		# Drop filename
		path = path[:path.rfind(os.path.sep)]
	elif path[-1] == os.path.sep:
		path = path[:-1]
	return path[path.rfind(os.path.sep)+1:]

def getDelim(line):
	# Returns delimiter from test file
	for i in ["\t", ",", " "]:
		if i in line:
			return i
	printError("Cannot determine delimeter. Check file formatting")

def __callProc__(cmd, sout, serr):
	# Calls process with given command and output
	try:
		call = Popen(split(cmd), stdout = sout, stderr = serr)
		call.wait()
		if call.returncode is not None:
			return True
		else:
			return False
	except:
		s = cmd.split()
		proc = s[0]
		if "-" not in s[1]:
			# Get subcommand if present
			proc += " " + s[1]
		elif proc == "java":
			# Replace call to jar with name of jar
			proc = getFileName(s[2])
		print(("\t[Warning] Could not call {}").format(proc), file=stderr)
		return False

def runProc(cmd, log = None):
	# Wraps call to Popen, writes stdout/stdout err to log/devnull/stdout&stderr, returns True if no errors
	if not log:
		log = os.devnull
	if log == "stdout" or log == "stderr":
		# Pipe output to stdout/err
		return __callProc__(cmd, stdout, stderr)
	else:
		with open(log, "w") as out:
			# Pipe output to log
			return __callProc__(cmd, out, out)
