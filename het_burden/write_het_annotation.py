"""
Categorize alleles in a hets/ folder and write output to hets/*/*.hets.annotations.

Looks for hets/cutoff/*.hets files in the give hets folder.
Outputs to the same location.
Example:
python categorize_alleles.py -d ./hets -l ./het_burden/lib

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


def read_hets_folder(hets_glob_str,
		gb_frequencies, mitobank, phylo, veps):
	for filename in glob.glob(hets_glob_str):
		output_filename = filename + '.annotations'
		read_hets_file_output_annotated_hets(filename, output_filename,
			gb_frequencies, mitobank, phylo, veps)


def read_hets_file_output_annotated_hets(filename, output_filename,
	gb_frequencies, mitobank, phylo, veps):
	hets_file = open(filename, 'r')
	annotations_file = open(output_filename, 'w')
	het_alleles = set()
	refseq = hettools.read_crs()
	header = "locus\tallele\tperson"
	for key in ['phylotree', 'genbank_freq', 'mitobank', 'mitobank_status',
			'sift', 'sift_status']:
		header += "\t%s" % key
	annotations_file.write("%s\n" % header)
	for li in hets_file:
		s = li.rstrip('\n').split('\t')
		locus = int(s[1])
		majAllele = s[4]
		majAlleleFreq = float(s[5])
		minAllele = s[6]
		minAlleleFreq = float(s[7])
		total_count = float(s[5]) + float(s[7])
		if total_count < 1:
			print "Error on line %s" % li
			continue
		fractionMaj = majAlleleFreq/total_count
		fractionMin = minAlleleFreq/total_count
		if(majAllele != refseq[locus-1]):
			allele_tuple = (locus, majAllele)
			annotations_maj = annotate_allele(allele_tuple,
				gb_frequencies, mitobank, phylo, veps)
			out_li = "%i\t%s\t%s%s\n" % (locus, majAllele, s[0],
				ann_line(annotations_maj))
			annotations_file.write(out_li)
		if(minAllele != refseq[locus-1]):
			allele_tuple = (locus, minAllele)
			annotations_min = annotate_allele(allele_tuple,
				gb_frequencies, mitobank, phylo, veps)
			out_li = "%i\t%s\t%s%s\n" % (locus, minAllele, s[0],
				ann_line(annotations_min))
			annotations_file.write(out_li)

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


def annotate_allele(allele_tuple,
		gb_frequencies, mitobank, phylo, veps):

	allele = allele_tuple
	not_in_mitobank = False
	annotations = {'phylotree': "Not in phylotree.",
			'genbank_freq': '0',
			'mitobank': "Not in mitobank.",
			'mitobank_status': "NA",
			'sift': "NA",
			'sift_status': 'NA'}
	if(allele in phylo):
		annotations['phylotree'] = str(phylo[allele])
	if allele in gb_frequencies:
		annotations['genbank_freq'] = str(gb_frequencies[allele])
	if(allele in mitobank):
		annotations['mitobank'] = mitobank[allele]
		if re.search('Reported', mitobank[allele]):
			annotations['mitobank_status'] = 'Reported'
		if re.search('Cfrm', mitobank[allele]):
			annotations['mitobank_status'] = 'Cfrm'
	if(allele in veps):
		annotations['sift'] = veps[allele]
		siftdel = re.search(r'SIFT=deleterious', siftline)
		if(siftdel is not None):
			annotations['sift_status'] = 'Deleterious.'
		else:
			annotations['sift_status'] = 'Not deletrious.'
	for key in annotations:
		annotations[key] = re.sub("\t", ",", annotations[key])
	return annotations

def process_vep_files(hets_glob_str):
	for filename in glob.glob(hets_glob_str):
		basename = os.path.basename(filename)
		print "hets filename:%s" % basename
		m = re.match("([a-zA-Z0-9]+)\.(.*)hets$", basename)
		if m is not None:
			keyword = m.group(1)
			vepfileguess = "%s.fullvep" % keyword
			print "guessing the vep file is %s" % vepfileguess
			if(os.path.exists(vepfileguess)):
				print "\tfound vep file"
				veps = hettools.readInterestingVepLines(vepfileguess)  # key by tuple
			else:
				print "\tNo such vep file"
				veps = {}
		else:
			print "\tno good guess for vep file"
			veps = {}
	return veps


if __name__ == '__main__':
	args = parse_args()
	in_folder = args.hets_folder
	for cutoff in ['high', 'mid', 'low']:
		print "\n***Cutoff: %s***" % cutoff
		hets_glob_str = in_folder + '/%s/*.hets' % cutoff
		gb_frequencies = hettools.getGenbankFrequencies()  # key by tuple
		mitobank = hettools.readMitobank()  # key by tuple
		phylo = hettools.read_phylotree()  # key by tuple
		observed_alleles = list()
		pathogenic_counts = dict()
		print "\n******************\n"
		if len(glob.glob(hets_glob_str)) == 0:
			print "No files match %s." % hets_glob_str
			sys.exit()
		veps = process_vep_files(hets_glob_str)
		read_hets_folder(hets_glob_str,
				gb_frequencies, mitobank, phylo, veps)

