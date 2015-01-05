import glob
import argparse
import os
import sys
import subprocess
import re
from operator import itemgetter


def parse_input():
	parser = argparse.ArgumentParser(description=
	"""Requires two arguments: an input folder of .bam files, 
	and an output folder for .bc files.
	For a given folder of .bam files, outputs .bc files
	containing the number of observations of each base from
	each direction at a given locus.
	""")
	parser.add_argument('-i', '--input_folder')
	parser.add_argument('-o', '--output_folder')
	parser.add_argument('-f', '--rcrs_file_path', 
		default='lib/rcrs.fa')
	args = parser.parse_args()
	# Check if input folder exists.
	if not os.path.exists(args.input_folder):
		print "No such input folder."
		sys.exit()
	# Check if output folder exists.
	if not os.path.exists(args.output_folder):
		print "No such input folder."
		sys.exit()
	return args


def parse_base(s):
	"""Parse a split line from an mpileup file 
	and return a line for a .bc file.
	mpileup file format is the following:
	chr loc ref depth calls quality
	This code is a rewrite of the Mitoseek Perl code in Python.
	"""
	calls = s[4]
	quals = s[5]
	calls = calls.translate(None, '$') #remove  $
	calls = re.sub( '\^.', '', calls)
	deletions = re.finditer( '-(\d+)', calls)
	if deletions is not None:
		for m in deletions:
			calls = re.sub( "-\d+\w{%i}" % int(m.group(1)), '', calls, count=1)
	insertions = re.finditer( '\+(\d+)', calls)
	if insertions is not None:
		for m in insertions:
			calls = re.sub( "\+\d+\w{%i}" % int(m.group(1)), '', calls, count=1)
	(fA, fT, fC, fG) = (0, 0, 0, 0)
	(rA, rT, rC, rG) = (0, 0, 0, 0)
	obsBases = { 'fA': 0, 'fT': 0, 'fC': 0, 'fG':0, 'rA': 0, 'rT': 0, 'rC': 0, 'rG': 0}
	for i in range(0, len(calls)): 
		try:
			if( ( ord(quals[i]) - 33 ) >= 20):
				if( calls[i] == 'A' ):
					obsBases['fA'] += 1 
				if( calls[i] == 'T' ):
					obsBases['fT'] += 1 
				if( calls[i] == 'C' ):
					obsBases['fC']+= 1 
				if( calls[i] == 'G' ):
					obsBases['fG'] += 1 
				if( calls[i] == 'a' ):
					obsBases['rA'] += 1 
				if( calls[i] == 't' ):
					obsBases['rT'] += 1 
				if( calls[i] == 'c' ):
					obsBases['rC'] += 1 
				if( calls[i] == 'g' ):
					obsBases['rG'] += 1
				if( calls[i] == '.' ):	
					if( s[2] == 'A' ):
						obsBases['fA'] += 1 
					if( s[2] == 'T' ):
						obsBases['fT'] += 1 
					if( s[2] == 'C' ):
						obsBases['fC'] += 1 
					if( s[2] == 'G' ):
						obsBases['fG'] += 1 
				if( calls[i] == ',' ):	
					if( s[2] == 'A' ):
						obsBases['rA'] += 1 
					if( s[2] == 'T' ):
						obsBases['rT'] += 1 
					if( s[2] == 'C' ):
						obsBases['rC'] += 1 
					if( s[2] == 'G' ):
						obsBases['rG'] += 1 
		except:
			print "Error parsing line %s." % str(s)
			print "Calls are %s with len %i." % (calls, len(calls))
			print "Quality scores are %s with len %i." % (
				quals, len(quals))
	l = "%s\t%s\t%s\t" % (s[0], str(s[1]), s[2])
	l += "%i\t%i\t%i\t%i\t" % (obsBases['fA'], obsBases['fT'], obsBases['fC'], obsBases['fG'])
	l += "%i\t%i\t%i\t%i\n" % (obsBases['rA'], obsBases['rT'], obsBases['rC'], obsBases['rG'])
	if(len(calls) == 0):
		l = "%s\t%s\t%s\t" % (s[0], str(s[1]), s[2])
		l += "%i\t%i\t%i\t%i\t" % (0, 0, 0, 0)
		l += "%i\t%i\t%i\t%i\n" % (0, 0, 0, 0)

	return l


def convert_all_bam_to_bc(bamfolder, bcfolder, rcrs_file_path):
	"""Converts every file in the input folder into a bc file.
	The bc file is written to the output folder.
	"""
	folderStr = bamfolder + '/*.bam'
	numfile = 0
	for name in glob.glob(folderStr):
		numfile += 1
		basename = re.search( r'/([_\wa-zA-Z\.0-9]+.bam)$', name).group(1)
		basecallfilename = bcfolder + basename + '.bc'
		if(os.path.exists(basecallfilename)):
			print "Basecall file %s exists; will not over-write." % (
				basecallfilename)
			continue  # Do not overwrite existing basecall files.
		outf = open(basecallfilename, 'w') 
		cmdl = [ "samtools", "mpileup", '-q', '20', '-f', 
		rcrs_file_path,  name ] 
		#  -q 20 means skip alignments with MAPQ<20. 
		#  This cutoff is the default for Mitosek.
		pileup = subprocess.Popen(cmdl, stdout=subprocess.PIPE)
		v = pileup.communicate()[0].split('\n')
		lastLoci = 0
		for a in v:
			if re.match(r'^#', a):
				continue
			s = a.rstrip('\n').split('\t')
			if len(s) == 6:
				if( int(s[1]) != (lastLoci + 1) ):
					for loci in range(lastLoci+1, int(s[1]), 1):
						l = "%s\t%s\t%s\t" % (s[0], loci, 'N')
						l += "%i\t%i\t%i\t%i\t" % (0, 0, 0, 0)
						l += "%i\t%i\t%i\t%i\n" % (0, 0, 0, 0)
						outf.write(l)
				lineO = parse_base(s) 
				outf.write(lineO)
				lastLoci = int(s[1])
			else:
				pass  # Line not expected length.
		outf.close()
	print "Processed %i files." % numfile


if __name__ == "__main__":
	args = parse_args()
	convert_all_bam_to_bc(args.input_folder, args.output_folder,
		args.rcrs_file_path)
