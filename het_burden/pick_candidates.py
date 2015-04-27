"""
Categorize alleles in a hets/ folder and write output to hets/*/*.hets.annotations.

Looks for hets/cutoff/*.hets files in the give hets folder.
Outputs to the same location.
Example:
python pick_candidates.py -d ./hets -l ./het_burden/lib

Requires pathogenicity information in the given lib/ folder (will
look in source directory for lib/).
"""
import re
import sys
import os
import hettools
import glob
import math
import argparse


def parse_args():
	src_path = os.path.dirname(os.path.realpath(__file__))
	parser = argparse.ArgumentParser(description=
	"""Process files indicating positions of heteroplasmy to generate concordance 
	and pathogenicity information.
	""")
	parser.add_argument('-d', '--hets_folder', default='none',
	help='(Required) hets_folder')
	parser.add_argument('-l', '--lib', default="%s/lib/" % src_path,
	help="""Library folder contianing map file, trios file, rcrs.fa, 
	mitobank_variant_freqs.txt, and whatever other information files 
	are required. The filenames given for these files will be expected
	in the lib folder.
	""")
	args = parser.parse_args() #args.justVcf now holds value
	print "Command line options given were:"
	return args


def is_rare(allele_tuple, gb_frequencies):
	allele = allele_tuple
	if(allele in gb_frequencies):
		observ = gb_frequencies[allele]
	else:
		observ = 0
	if(observ < 6):
		return True
	else:
		return False


def read_hets_folder(glob_str):
	if len(glob.glob(glob_str)) == 0:
		print "No files match %s." % glob_str
		sys.exit()
	for filename in glob.glob(glob_str):
		if re.search('filled_in', filename) is not None:
			continue
		annotated_filename = filename + '.annotations'
		print "Reading .hets file %s..." % filename
		hets = read_hets_file(filename)
		print "Reading .annotated file %s..." % annotated_filename
		ann_hets = read_annotated_hets_file(annotated_filename, hets)
		print "Prioritizing..."
		candidates = find_serious_alleles(ann_hets)
		candidates_filename = filename + '.candidates'
		print "Writing %s..." % candidates_filename
		print_candidates(candidates, candidates_filename)


def read_hets_file(filename):
	hets = {}
	hets_file = open(filename, 'r')
	for li in hets_file:
		s = li.rstrip('\n').split('\t')
		person = s[0]
		locus = int(s[1])
		hets.setdefault(person, {})
		hets[person][locus] = li
#		majAllele = s[4]
#		majAlleleFreq = float(s[5])
#		minAllele = s[6]
#		minAlleleFreq = float(s[7])
#		fractionMaj = majAlleleFreq/total_count
#		fractionMin = minAlleleFreq/total_count
		total_count = float(s[5]) + float(s[7])
		if total_count < 1:
			print "Error on line %s" % li
			continue
	return hets


def read_annotated_hets_file(filename, hets):
	hets_file = open(filename, 'r')
#	refseq = hettools.read_crs()
#	header = "locus\tallele\tperson"
#	for key in ['phylotree', 'genbank_freq', 'mitobank', 'mitobank_status',
#			'sift', 'sift_status']:
#		header += "\t%s" % key
	header = next(hets_file).rstrip('\n').split('\t')
	annotations = header[3:len(header)]
	ann_hets = {}
	for ann_li in hets_file:
		s = ann_li.rstrip('\n').split('\t')
		person = s[2]
		if person not in hets:
			print "Error: unexpected person %s in .annotated file." % person
		locus = int(s[0])
		if locus not in hets[person]:
			print "Error: unexpected locus %i on line %s" % (locus, ann_li)
		allele = s[1]
		allele_info = s[3:len(s)]
		allele_dict = dict(zip(annotations, allele_info))
		ann_hets.setdefault(person, {})
		ann_hets[person][locus] = {'hets_line': hets[person][locus],
					'ann_line': ann_li,
					'allele_info': allele_dict}
		#print allele_dict
	return ann_hets


def print_candidates(candidates, filename):
	print "print_candidates(), filename=%s" % filename
	with open(filename, 'w') as f:
		for person in candidates:
			for locus in candidates[person]:
	#			li = candidates[person][locus]['hets_line'].rstrip('\n')
				li = candidates[person][locus]['ann_line'].rstrip('\n')
				li += candidates[person][locus]['hets_line']
	#				candidates[person][locus]['allele_info'] + "\n")
				f.write(li)
				print "Printing the following line:"
				print li.rstrip('\n')

def find_serious_alleles(ann_hets):
	candidates = {}
	for person in ann_hets:
		for locus in ann_hets[person]:
			if is_serious(ann_hets[person][locus]):
				candidates.setdefault(person, {})
				candidates[person][locus] = ann_hets[person][locus]
	return candidates

def is_serious(allele_dict):
	allele_info = allele_dict['allele_info']
	hets_line = allele_dict['hets_line']
	try:
		frequency = float(hets_line.split('\t')[2])
	except:
		print "Could not get frequency from 'hets_line' in %s" % str(allele_dict)
		print "Specifically the value: %s" % str(allele_dict['hets_line'])
		frequency = 0
		sys.exit()
	if frequency < 0.3:
		return False
	is_patho = False
	if allele_info['mitobank_status'] == 'Cfrm':
		is_patho = True
		return True
	return False


def ann_line(annotations):
	"""Create a string from annotations dict for writing to annotated hets file.
	"""
	li = ""
	for key in ['phylotree', 'genbank_freq', 'mitobank', 'mitobank_status',
			'sift', 'sift_status']:
		if key not in annotations:
			print "Did not find %s in %s..." % (key, str(annotations))
			sys.exit()
		li += "\t%s" % str(annotations[key])
	return li


if __name__ == '__main__':
	args = parse_args()
	in_folder = args.hets_folder
	for cutoff in ['high', 'mid', 'low']:
		print "\n***Pick candidates, cutoff: %s***" % cutoff
		hets_glob_str = in_folder + '/%s/*.hets' % cutoff
		read_hets_folder(hets_glob_str)
		hets_glob_str = in_folder + '/%s/*.homo' % cutoff
		read_hets_folder(hets_glob_str)
		hets_glob_str = in_folder + '/%s/*.other' % cutoff
		read_hets_folder(hets_glob_str)

