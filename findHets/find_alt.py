"""
Finds heteroplasmic SNPs, homoplasmic SNPs, and SNPs that are majority non-reference
but do not fall into the categories of heteroplasmic or homoplasmic.

Example:
	python find_alt.py -b 1kgenomesbc -o 1kgenomes

Outputs .hets, .homo, and .other files in ./hets/cutoff.
Also outputs .vcf and .stats files for heteroplasmies.

"""
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
	"""Finds heteroplasmic and homoplasmic SNPs from a folder of
	basecall files.
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


def find_alt_loci_in_bc_files(folder, refseq):
	"""Identify heteroplasmic loci in every basecall file in a folder.
	folder -- string in format /path/*bc
	Checks if loci passes filters, then sets the het dict if it does.
	het key=basename of .bam file. value = dict()
	het[basename] key = locus. value = sorted list of base observations tuple.
	"""
	min_total_reads_by_strand = 5
	min_count_allele_by_strand = 1
	het = dict()
	homoplasmic = dict()
	majority_alt = dict()
	baseToIndexF = {'A': 3, 'T': 4, 'C': 5, 'G': 6}
	baseToIndexR = {'A': 7, 'T': 8, 'C': 9, 'G': 10}
	for fname in glob.glob(folder):
		basename = os.path.basename(fname).rstrip('.bc')
		het[basename] = dict()
		homoplasmic[basename] = dict()
		majority_alt[basename] = dict()
		bcfile = open(fname, 'r')
		for li in bcfile:
			q = li.rstrip('\n').split('\t')
			s = q[0:3]
			for i in q[3:11]:
				s.append(int(i))
			totalfor = sum(s[3:7])
			totalrev = sum(s[7:11])
			if((totalfor < min_total_reads_by_strand)
				or (totalrev < min_total_reads_by_strand)):
				# This site not is usable.	
				continue
			obsAlleles = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
			for i, a in enumerate(['A', 'T', 'C', 'G'], start=3):
				obsAlleles[a] = s[i] + s[4+i]
			sorted_bc = sorted(obsAlleles.items(), key=itemgetter(1))
			majAllele = sorted_bc[-1][0]	
			minorAllele = sorted_bc[-2][0]
			if majAllele != refseq[int(s[1])-1] and (
				sorted_bc[-2][1] == 0):
				homoplasmic[basename][s[1]] = li
			if majAllele != refseq[int(s[1])-1] and (
				sorted_bc[-2]):
				majority_alt[basename][s[1]] = li
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
	print """Het loci: %i, Homo loci: %i, Majority alt loci: %i""" % (
		sum([len(het[person]) for person in het]),
		sum([len(homoplasmic[person]) for person in homoplasmic]),
		sum([len(majority_alt[person]) for person in majority_alt]))
	return (het, homoplasmic, majority_alt)


def filter_alts(cutoff, alts):
	"""Write hets file indicating where heteroplasmic loci are,
	as well as .vcf files.
	Uses the het dict set in find_het_loci_in_bc_files. 
	"""
	baseToIndexF = {'A': 3, 'T': 4, 'C': 5, 'G': 6}
	baseToIndexR = {'A': 7, 'T': 8, 'C': 9, 'G': 10}
	if(cutoff=='high'):
		min_total_reads_by_strand = 10
	if(cutoff=='mid'):
		# mid has a special OR condition.
		min_total_reads_high = 10
		min_total_reads_low = 5
	if(cutoff=='low'):
		min_total_reads_by_strand = 5
	filtered = dict()
	for person in alts:
		for locus in alts[person]:
			res_list = line_to_sorted_bc(alts[person][locus])
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
			if(cutoff=='high' or cutoff=='low'):
				if (not((totalfor>=min_total_reads_by_strand) 
					and (totalrev>=min_total_reads_by_strand))):
					continue 
			filtered.setdefault(person, dict())
			filtered[person][locus] = sorted_bc
	return filtered


def filter_and_write_alts(cutoff, hets, homoplasmic, majority_alt, basename, refseq):
	if(not os.path.exists(r'./hets')):
		os.system(r'mkdir hets')
	if(not os.path.exists(r'./hets/high')):
		os.system(r'mkdir hets/high')
	if(not os.path.exists(r'./hets/mid')):
		os.system(r'mkdir hets/mid')
	if(not os.path.exists(r'./hets/low')):
		os.system(r'mkdir hets/low')
	(filtered_hets, num_by_maf) = write_hets(cutoff, hets, out_basename, refseq)
	if(cutoff == 'high'):
		outF = open("hets/%s/%s.homo" % (cutoff, out_basename), 'w')
	else:
		outF = open("hets/%s/%s.homo" % (cutoff, out_basename), 'w')
	filtered_homoplasmic = filter_alts(cutoff, homoplasmic)
	filtered_homoplasmic = subtract_loci(filtered_homoplasmic, filtered_hets)
	write_alts(filtered_homoplasmic, outF, write_freq_as_maj_allele=True)
	if(cutoff == 'high'):
		outF = open("hets/%s/%s.other" % (cutoff, out_basename), 'w')
	else:
		outF = open("hets/%s/%s.other" % (cutoff, out_basename), 'w')
	filtered_majority_alt = filter_alts(cutoff, majority_alt)
	filtered_majority_alt = subtract_loci(filtered_majority_alt, filtered_hets)
	filtered_majority_alt = subtract_loci(filtered_majority_alt, filtered_homoplasmic)
	write_alts(filtered_majority_alt, outF, write_freq_as_maj_allele=True)
	return (filtered_hets, num_by_maf)


def subtract_loci(dict1, dict2):
	filtered_dict1 = dict()
	for person in dict1:
		if person not in dict2:
			filtered_dict1[person] = dict1[person]
			continue
		for loci in dict1[person]:
			if loci not in dict2[person]:
				filtered_dict1.setdefault(person, {})
				filtered_dict1[person][loci] = dict1[person][loci]
				continue
	return filtered_dict1


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
	outF = open("hets/%s/%s.hets" % (cutoff, out_basename), 'w')
	num_by_maf = write_alts(filtered_het, outF)
	return (filtered_het, num_by_maf)


def write_alts(filtered_loci, outF, write_freq_as_maj_allele=False):
	filtered_het = filtered_loci
	num_by_maf = dict()
	for index in range(0, 60, 1):
		fraction = float(index)/100.0
		num_by_maf[fraction] = 0
	for person in filtered_het:
		largestFraction = 0
		for locus in filtered_het[person]:
			depth = sum([x[1] for x in filtered_het[person][locus]])
			fractionMinor = float(filtered_het[person][locus][-2][1])/float(depth)
			if(largestFraction < fractionMinor):
				largestFraction = fractionMinor
			if write_freq_as_maj_allele:
				frequency = float(
					filtered_het[person][locus][-1][1])/float(depth)
			else:
				frequency = largestFraction
			outF.write("%s\t%s\t%f\t%i\t%s\t%f\t%s\t%f\t%s" % (
				person, str(locus), frequency, depth,
				filtered_het[person][locus][-1][0],
				float(filtered_het[person][locus][-1][1]), 
				filtered_het[person][locus][-2][0],
				float(filtered_het[person][locus][-2][1]),
				str(filtered_het[person][locus])))
			if(hettools.is_bad_loci(locus)):
				outF.write("\tLoci is at bad allele.")
			outF.write("\n")
		largestFraction = float(math.trunc(largestFraction * 100))/float(100)
		num_by_maf[largestFraction] += 1
	outF.close()
	return (num_by_maf)


def write_vcf(het, cutoff, refseq, no_keywords=False):
	if(not no_keywords):
		allVcf = open("hets/%s/%s.vcf" % (cutoff, out_basename), 'w')
		majorVcf = open("hets/%s/%s.hetMajor.vcf" % (cutoff, out_basename), 'w')
		minorVcf = open("hets/%s/%s.hetMinor.vcf" % (cutoff, out_basename), 'w')
	if(no_keywords):
		allVcf = open("hets/%s/all.vcf" % (cutoff, cutoff), 'w')
		majorVcf = open("hets/%s/all.hetMajor.vcf" % (cutoff, cutoff), 'w')
		minorVcf = open("hets/%s/all.hetMinor.vcf" % (cutoff, cutoff), 'w')
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
	outF = open("hets/%s/%s.hets.stats" % (cutoff, out_basename), 'w')
	lineO = "Number of individuals: %i" % len(het)
	for frequency in sorted(num_by_maf.keys()):
		lineO += "\nFraction with a highest abundance heteroplasmy between "
		if len(het) > 0:
			ratio = float(num_by_maf[frequency])/float(len(het))
		else:
			ratio = 0
		lineO += " %i%% and %i%%:" % (int(frequency*100), int(frequency*100+1))
		lineO += " %f" % ratio
	outF.write(lineO)
	outF.close()


def apply_cutoff(cutoff, hets, homoplasmic, majority_alt, out_basename, refseq):
	#(het, num_by_maf) = write_hets(cutoff, het, out_basename, refseq)
	(filtered_hets, num_by_maf) = filter_and_write_alts(
		cutoff, hets, homoplasmic, majority_alt, out_basename, refseq)
	write_vcf(filtered_hets, cutoff, refseq)
	write_stats(cutoff, num_by_maf, filtered_hets, out_basename)


if __name__ == '__main__':
	args = parse_input()
	src_path = os.path.dirname(os.path.realpath(__file__))
	if(not args.output):
		out_basename = "default_prefix"
	else:
		out_basename = args.output
	refseq = hettools.read_crs("%s/lib/rcrs.fa" % src_path)
	(het, homoplasmic, majority_alt) = find_alt_loci_in_bc_files(args.bc_folder + '*.bc', refseq)
	apply_cutoff('high', het, homoplasmic, majority_alt, out_basename, refseq)
	apply_cutoff('mid', het, homoplasmic, majority_alt, out_basename, refseq)
	apply_cutoff('low', het, homoplasmic, majority_alt, out_basename, refseq)
