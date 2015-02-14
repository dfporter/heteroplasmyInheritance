import argparse
import glob
import os
import re
from operator import itemgetter
import sys
import math
import hettools


def parse_input():
	parser = argparse.ArgumentParser(description=
	"""Process files indicating positions of heteroplasmy to generate concordance 
	and pathogenicity information.
	""")
	parser.add_argument('-b', '--bc_folder', 
	help='Folder of basecall files to parse.')
	parser.add_argument('-o', '--output', default='default_prefix',
	help='Optional prefix for output files.')
	#parser.add_argument('-p', '--paths',
	#help='Paths file.', default='./lib/paths2.paths')
	args = parser.parse_args()
	args.bc_folder = args.bc_folder.rstrip('/') + '/'
	args.output = args.output.rstrip('/') 
	if(not os.path.exists(args.bc_folder)):
		print ".bc folder does not exist"
	return args


def line_to_sorted_bc(li):
	"""Return an array and sorted list of base observations.
	Argument is a line from a basecall file. Returns a list of 
	[(array representation of line, with numbers converted to integers),
	(sorted list of tuples in the form ('A', 23), representing observations
	of alleles].
	"""
	totaldepth = 0
	totalfor = 0
	totalrev = 0
	q = li.rstrip('\n').split('\t')
	s = q[0:3]
	for i in q[3:11]:
		s.append(int(i))
	totalfor = sum(s[3:7])
	totalrev = sum(s[7:11])
	obsAlleles = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
	for i, a in enumerate(['A', 'T', 'C', 'G'], start=3):
		obsAlleles[a] = s[i] + s[4+i]
	sorted_bc = sorted(obsAlleles.items(), key=itemgetter(1))
	return [s, sorted_bc]


def find_het_loci_in_bc_files(folder):
	"""Identify heteroplasmic loci in every basecall file in a folder.
	folder -- string in format /path/*bc
	Checks if loci passes filters, then sets the het dict if it does.
	het key=basename of .bam file. value = dict()
	het[basename] key = locus. value = sorted list of base observations tuple.
	"""
	min_total_reads_by_strand = 5
	min_count_allele_by_strand = 1
	het = dict()
	baseToIndexF = {'A': 3, 'T': 4, 'C': 5, 'G': 6}
	baseToIndexR = {'A': 7, 'T': 8, 'C': 9, 'G': 10}
	for fname in glob.glob(folder):
		basename = os.path.basename(fname).rstrip('.bc')
		het[basename] = dict()
		bcfile = open(fname, 'r')
		for li in bcfile:
			totaldepth = 0
			totalfor = 0
			totalrev = 0
			q = li.rstrip('\n').split('\t')
			s = q[0:3]
			for i in q[3:11]:
				s.append(int(i))
			totalfor = sum(s[3:7])
			totalrev = sum(s[7:11])
			obsBases = dict()
			if((totalfor >= min_total_reads_by_strand)
				and (totalrev >= min_total_reads_by_strand)): 
				# This site is usable.	
				obsAlleles = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
				for i, a in enumerate(['A', 'T', 'C', 'G'], start=3):
					obsAlleles[a] = s[i] + s[4+i]
				sorted_bc = sorted(obsAlleles.items(), key=itemgetter(1))
				majAllele = sorted_bc[-1][0]	
				minorAllele = sorted_bc[-2][0]	
				if((s[baseToIndexF[majAllele]] >= min_count_allele_by_strand) 
				and (s[baseToIndexR[majAllele]] >= min_count_allele_by_strand) 
				and (s[baseToIndexF[minorAllele]] >= min_count_allele_by_strand) 
				and (s[baseToIndexR[minorAllele]] >= min_count_allele_by_strand)):
					# We are heteroplasmic
					het[basename][s[1]] = li
				else:
					# No heteroplasmy
					pass 
		bcfile.close()
	return het


def write_hets(cutoff, het, out_basename, refseq):
	"""Write hets file indicating where heteroplasmic loci are,
	as well as .vcf files.
	Uses the het dict set in find_het_loci_in_bc_files. 
	"""
	baseToIndexF = {'A': 3, 'T': 4, 'C': 5, 'G': 6}
	baseToIndexR = {'A': 7, 'T': 8, 'C': 9, 'G': 10}
	if(not os.path.exists(r'./hets')):
		os.system(r'mkdir hets')
	if(not os.path.exists(r'./hets/high')):
		os.system(r'mkdir hets/high')
	if(not os.path.exists(r'./hets/mid')):
		os.system(r'mkdir hets/mid')
	if(not os.path.exists(r'./hets/low')):
		os.system(r'mkdir hets/low')
	if(cutoff=='high'):
		min_total_reads_by_strand = 10
		min_allele_count_by_strand = 2
	if(cutoff=='mid'):
		# mid has a special OR condition.
		min_total_reads_high = 10
		min_total_reads_low = 5
		min_count_of_allele_high = 2
		min_count_of_allele_low = 1
	if(cutoff=='low'):
		min_total_reads_by_strand = 5
		min_allele_count_by_strand = 1
	filtered_het = dict()
	for p in het:
		for locus in het[p]:
			res_list = line_to_sorted_bc(het[p][locus])
			s = res_list[0]
			sorted_bc = res_list[1]
			totalfor = sum(s[3:7])
			totalrev = sum(s[7:11])
			majAllele = sorted_bc[-1][0]	
			minorAllele = sorted_bc[-2][0]
			if(cutoff=='mid'):
				if (not((totalfor>=min_total_reads_low) 
					and (totalrev>=min_total_reads_low))):
					continue 
				if (not((totalfor>=min_total_reads_high) 
					or (totalrev>=min_total_reads_high))):
					continue 
				if( (s[baseToIndexF[majAllele]] >= min_count_of_allele_low) 
				and (s[baseToIndexR[majAllele]] >= min_count_of_allele_low) 
				and (s[baseToIndexF[minorAllele]] >= min_count_of_allele_low) 
				and (s[baseToIndexR[minorAllele]] >= min_count_of_allele_low)
				# The above block ensures both major and minor alleles are above
				# the lower cutoff.
				and (s[baseToIndexF[majAllele]] >= min_count_of_allele_high
				or s[baseToIndexR[majAllele]] >= min_count_of_allele_high)
				and (s[baseToIndexF[minorAllele]] >= min_count_of_allele_high
				or s[baseToIndexR[minorAllele]] >= min_count_of_allele_high)
				# The above block ensures that at least one strand direction
				# passes the higher cutoff.
					):
					# Locus is heteroplasmic.
					pass
				else:
					continue
			if(cutoff=='high' or cutoff=='low'):
				if (not((totalfor>=min_total_reads_by_strand) 
					and (totalrev>=min_total_reads_by_strand))):
					continue 
				if( (s[baseToIndexF[majAllele]] >= min_allele_count_by_strand) 
				and (s[baseToIndexR[majAllele]] >= min_allele_count_by_strand) 
				and (s[baseToIndexF[minorAllele]] >= min_allele_count_by_strand) 
				and (s[baseToIndexR[minorAllele]] >= min_allele_count_by_strand)):
				# Locus is heteroplasmic.
					pass
				else:
					continue
			filtered_het.setdefault(p, dict())
			filtered_het[p][locus] = sorted_bc
	het = filtered_het
	numberOfIndividualsAboveMaf = {0: 0, 0.1: 0, 0.2: 0, 0.3: 0, 0.4: 0, 0.5: 0, 0.6: 0}
	if(cutoff == 'high'):
		outF = open("hets/%s/%s.hets" % (cutoff, out_basename), 'w')
	else:
		outF = open("hets/%s/%s.%s.hets" % (cutoff, out_basename, cutoff), 'w')
	for person in het:
		largestFraction = 0
		for locus in het[person]:
			depth = sum([x[1] for x in het[person][locus]])
			fractionMinor = float(het[person][locus][-2][1])/float(depth)
			if(largestFraction < fractionMinor):
				largestFraction = fractionMinor
			outF.write("%s\t%s\t%f\t%i\t%s\t%f\t%s\t%f\t%s" % (
				person, str(locus), fractionMinor, depth, 
				het[person][locus][-1][0], float(het[person][locus][-1][1]), 
				het[person][locus][-2][0], float(het[person][locus][-2][1]), 
				str(het[person][locus])))
			if(hettools.is_bad_loci(locus)):
				outF.write("\tLoci is at bad allele.")
			outF.write("\n")
		largestFraction = float(math.trunc(largestFraction * 10))/float(10)
		numberOfIndividualsAboveMaf[largestFraction] += 1
	outF.close()
	return (het, numberOfIndividualsAboveMaf)


def write_vcf(het, cutoff, refseq, no_keywords=False):
	if(cutoff == 'high' and not no_keywords):
		allVcf = open("hets/%s/%s.vcf" % (cutoff, out_basename), 'w')
		majorVcf = open("hets/%s/%s.hetMajor.vcf" % (cutoff, out_basename), 'w')
		minorVcf = open("hets/%s/%s.hetMinor.vcf" % (cutoff, out_basename), 'w')
	if(cutoff != 'high' and not no_keywords):
		allVcf = open("hets/%s/%s.%s.vcf" % (cutoff, out_basename, cutoff), 'w')
		majorVcf = open("hets/%s/%s.%s.hetMajor.vcf" % (cutoff, out_basename, cutoff), 'w')
		minorVcf = open("hets/%s/%s.%s.hetMinor.vcf" % (cutoff, out_basename, cutoff), 'w')
	if(no_keywords):
		allVcf = open("hets/%s/all.%s.vcf" % (cutoff, cutoff), 'w')
		majorVcf = open("hets/%s/all.%s.hetMajor.vcf" % (cutoff, cutoff), 'w')
		minorVcf = open("hets/%s/all.%s.hetMinor.vcf" % (cutoff, cutoff), 'w')
	lines_written = set()
	for p in het:
		for locus in het[p]:
			try:
				depth = sum([x[1] for x in het[p][locus]])
			except:
				print "Problem with p=%s, locus=%s: %s" % (str(p), str(locus),
					str(het[p][locus]))
				continue
			fractionMinor = float(het[p][locus][-2][1])/float(depth)
 			# Write major allele vcf.
			if(het[p][locus][-1][0] != refseq[int(locus)-1]):
				vcfLine = "MT\t%i\t.\t%s\t%s\t999\t.\tDP=1000\n" % (
				int(locus), refseq[int(locus)-1] , het[p][locus][-1][0])
				majorVcf.write(vcfLine)
				if(vcfLine not in lines_written):
					allVcf.write(vcfLine)
					lines_written.add(vcfLine)
			# Write minor allele vcf.
			if(het[p][locus][-2][0] != refseq[int(locus)-1]):
				vcfLine = "MT\t%i\t.\t%s\t%s\t999\t.\tDP=1000\n" % (
				int(locus), refseq[int(locus)-1], het[p][locus][-2][0]) 
				minorVcf.write(vcfLine)
				if(vcfLine not in lines_written):
					allVcf.write(vcfLine)
					lines_written.add(vcfLine)
	allVcf.close()
	majorVcf.close()
	minorVcf.close()


def write_stats(cutoff, num_by_maf, het, out_basename):
	"""Writes the number of individuals above each maf cutoff"""
	if(cutoff == 'high'):
		outF = open("hets/%s/%s.hets.stats" % (cutoff, out_basename), 'w')
	else:
		outF = open("hets/%s/%s.%s.hets.stats" % (cutoff, out_basename, cutoff), 'w')
	lineO = "Number of individuals: %i" % len(het)
	for frequency in sorted(num_by_maf.keys()):
		lineO += "\nFraction with a highest abundance heteroplasmy between "
		if len(het) > 0:
			ratio = float(num_by_maf[frequency])/float(len(het))
		else:
			ratio = 0
		lineO += " %i%% and %i%%:" % (int(frequency*100), int(frequency*100+10))
		lineO += " %f" % ratio
	outF.write(lineO)
	outF.close()


def apply_cutoff(cutoff, het, out_basename, refseq):
	(het, num_by_maf) = write_hets(cutoff, het, out_basename, refseq)
	write_vcf(het, cutoff, refseq)
	write_stats(cutoff, num_by_maf, het, out_basename)


if __name__ == '__main__':
	args = parse_input()
	src_path = os.path.dirname(os.path.realpath(sys.argv[0]))
	if(not args.output):
		out_basename = "default_prefix"
	else:
		out_basename = args.output
	refseq = hettools.read_crs("%s/lib/rcrs.fa" % src_path)
	het = find_het_loci_in_bc_files(args.bc_folder + '*.bc')
	apply_cutoff('high', het, out_basename, refseq)
	apply_cutoff('mid', het, out_basename, refseq)
	apply_cutoff('low', het, out_basename, refseq)
