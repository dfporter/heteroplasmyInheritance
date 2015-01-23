"""
Categorize alleles in a hets/ folder and write output to terminal.

Example:
python categorize_alleles.py -d ./hets -l ./het_burden/lib

No files are written.

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
NUM_IND = 12577


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


def print_obs_alleles(observed_alleles):
	obs = observed_alleles
	bins = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
	# values above 4 should be 0
	highest_in_ind = dict()
	# list of (file, alleletuple, fraction )
	for anobs in obs:
		anobs[2] = math.floor(10*anobs[2])
		if(anobs[2] == 5):
			print "Error: minor allele set 0.5-.6."
			print anobs
			anobs[2] = 4
		# 7% = 0. 23% = 2. 45%=4.
		if anobs[0] in highest_in_ind:
			if(anobs[2] > highest_in_ind[anobs[0]]):
				highest_in_ind[anobs[0]] = anobs[2]
		else:
			highest_in_ind[anobs[0]] = anobs[2]
	for anind in highest_in_ind:
		bins[highest_in_ind[anind]] += 1
	abovebin = dict()
	for abin in bins:
		abovebin[0] = bins[0] + bins[1] + bins[2] + bins[3] + bins[4] + bins[5]
		abovebin[1] = bins[1] + bins[2] + bins[3] + bins[4] + bins[5]
		abovebin[2] = bins[2] + bins[3] + bins[4] + bins[5]
		abovebin[3] = bins[3] + bins[4] + bins[5]
		abovebin[4] = bins[4] + bins[5]
		abovebin[5] = bins[5]
	li = "% het.\tn above given %"
	li += "\t% above given %het"
	li += "\tn in given bin"
	li += "\t% in given bin"
	# non-het individuals:
	li += "\nnot_het\t%i\t%.0f" % (NUM_IND, 100.0 * float(1))
	li += "\t%i" % (NUM_IND-len(highest_in_ind))
	li += "\t%.2f" % (100.0 * float(NUM_IND-len(highest_in_ind))/float(NUM_IND))
	for abin in bins:
		li += "\n%.0f" % float(10*abin)
		li += "\t%i\t%.2f" % (abovebin[abin],
					100.0 * float(abovebin[abin])/float(NUM_IND))
		li += "\t%i\t%.2f" % (bins[abin], 100.0 * float(bins[abin])/float(NUM_IND))
	print li


def print_pathogenic_counts(cat, mitobank, veps,
	pathogenic_counts, alt_format=False, write_ind=True):
	line = "\nAlleles listed in mitobank:"
	for allele in pathogenic_counts:
		mitob = str(mitobank[allele]).split('\t')
		if(alt_format):
			line += "\n%s & %s & %s & %i \\\\" % (
				mitob[3], mitob[2], mitob[8],
				len(pathogenic_counts[allele]))
		else:
			mitob = "%s\t%s" % ("\t".join(mitob[0:4]), mitob[8])
			mitob = str(mitob)
			line += "\n%s" % mitob
			line += "\tindividuals:%i" % len(pathogenic_counts[allele])
		if (not write_ind):
			continue
		for ind in pathogenic_counts[allele]:
			line += "\t%s:%s" % (ind[0], str(ind[1]))
	print line
	print "end print pathogenic"


def print_cat(cat, mitobank, veps):
	dumb = '''
	cat['n_all'] = len(list(cat['all']))
	cat['n_mitobank'] = len(list(cat['mitobank']))
	cat['n_mitobank_common'] = len(list(cat['mitobank_common']))
	cat['n_mitobank_rare'] = len(list(cat['mitobank_rare']))
	cat['n_phylotree'] = len(list(cat['phylotree']))
	cat['n_novel'] = len(list(cat['novel']))
	cat['n_novel_deleterious'] = len(list(cat['novel_deleterious']))
	cat['n_novel_nondeleterious'] = len(list(cat['novel_nondeleterious']))
	cat['n_not_mitobank_but_observed'] = len(list(cat['not_mitobank_but_observed']))
	cat['n_not_mitobank_but_observed_and_rare'] = len(list(cat['not_mitobank_but_observed_and_rare']))
	cat['n_not_mitobank_but_observed_and_common'] = len(list(cat['not_mitobank_but_observed_and_common']))
	cat['n_not_mitobank_but_observed_deleterious'] = len(list(cat['not_mitobank_but_observed_deleterious']))
	cat['n_not_mitobank_but_observed_nondeleterious'] = len(list(cat['not_mitobank_but_observed_nondeleterious']))
	cat['n_homoplasmy_change'] = len(list(cat['homoplasmy_change']))
'''
	cat_with_n = {}
	for key in cat:
		key_for_n = "n_" + key
		cat_with_n[key_for_n] = len(list(cat[key]))
	cat.update(cat_with_n)
	# build output
	line = "All variant alleles considered: %i" % cat['n_all']
	line += "\nAll variant alleles at a minimum of 10 reads: %i" % cat['n_all_above_depth']
	line += "\nCharacteristic of a reference haplotype: %i" % cat['n_phylotree']
	line += "\nConsidered deleterious (in mitobank): %i" % cat['n_mitobank']
	line += "\n	* Common: %i " % cat['n_mitobank_common']
	line += "\n	* Rare: %i " % cat['n_mitobank_rare']
	line += "\nPreviously observed, and not annotated as deleterious by mitobank:"
	line += " %s" % cat['n_not_mitobank_but_observed']
	line += "\n 	* Not predicted to be deleterious:"
	line += " %i" % cat['n_not_mitobank_but_observed_nondeleterious']
	line += "\n	* Predicted to be deleterious:"
	line += " %i" % cat['n_not_mitobank_but_observed_deleterious']
	line += "\n	* Common: %i" % cat['n_not_mitobank_but_observed_and_common']
	line += "\n	* Rare: %i" % cat['n_not_mitobank_but_observed_and_rare']
	line += "\nNot in GenBank: %i" % cat['n_novel']
	line += "\n	* Not predicted to be deleterious: %i" % cat['n_novel_nondeleterious']
	line += "\n 	* Predicted to be deleterious: %i" % cat['n_novel_deleterious']
	line += "\nDiscordance: "
	line += "\n	* Almost a switch in homoplasmy (max shared 0.2):"
	line += " %i" % cat['n_almost_homoplasmy_change']
	line += "\n	* True change in homoplasmy: %i " % cat['n_homoplasmy_change']
	line += "\n"
	print line
	line = "\nRare alleles listed in mitobank:"
	for allele in cat['mitobank_rare']:
		line += "\n%s" % str(mitobank[allele])
	line += "\nCommon alleles listed in mitobank:"
	for allele in cat['mitobank_common']:
		line += "\n%s" % str(mitobank[allele])
	line += "\nPreviously observed, not deleterious in mitobank, and:"
	line += "\n\tnot predicted to be deleterious:"
	line += "\t%i" % cat['n_not_mitobank_but_observed_nondeleterious']
	line += "\n\tpredicted to be deleterious:"
	for allele in cat['not_mitobank_but_observed_deleterious']:
		if(allele in veps):
			line += "\n%s:%s" % (str(allele), str(veps[allele]))
		else:
			line += "\n%s" % str(allele)
	line += "\nNot previously observed (in genbank) and:"
	line += "\n\tnot predicted to be deleterious:"
	line += "\t%i" % cat['n_novel_nondeleterious']
	line += "\n\tpredicted to be deleterious:"
	for allele in cat['novel_deleterious']:
		if(allele in veps):
			line += "\n%s:%s" % (str(allele), str(veps[allele]))
		else:
			line += "\n%s" % str(allele)
	line += "\n"
	print line


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
	gb_frequencies, mitobank, phylo, veps, cat,
	observed_alleles=[], pathogenic_counts=dict()):
	for filename in glob.glob(hets_glob_str):
		read_hets_format(filename,
			gb_frequencies, mitobank, phylo, veps, cat,
			observed_alleles=[], pathogenic_counts=dict())


def read_hets_format(filename,
	gb_frequencies, mitobank, phylo, veps, cat,
	observed_alleles=[], pathogenic_counts=dict()):
	f = open(filename, 'r')
	het_alleles = set()
	refseq = hettools.read_crs()
	for li in f:
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
		observed_alleles.append([s[0], (locus, minAllele), fractionMin])
		if(majAllele != refseq[locus-1]):
			allele_tuple = (locus, majAllele)
			annotate_allele(allele_tuple,
				gb_frequencies, mitobank, phylo, veps, cat,
				pathogenic_counts=pathogenic_counts,
				fraction=majAlleleFreq/total_count,
				individual=s[0])
		if(minAllele != refseq[locus-1]):
			allele_tuple = (locus, minAllele)
			annotate_allele(allele_tuple,
				gb_frequencies, mitobank, phylo, veps, cat,
				pathogenic_counts=pathogenic_counts,
				fraction=minAlleleFreq/total_count,
				individual=s[0])


def read_concordance_checked_format(filename_in,
	gb_frequencies, mitobank, phylo, veps, cat, observed_alleles=[]):
	"""
	This format is:
	trio_name	313G	max_change	...	mbc	A	3	ect.
	"""
	f = open(filename_in, 'r')
	het_alleles = set()
	for li in f:
		s = li.rstrip('\n').split('\t')
		if(len(s) < 3):
			continue
		cat['all'] += 1
		m = re.match(r'\A(\d+)([a-zA-Z]+)', s[1])
		locus = int(m.group(1))
		mut = m.group(2)
		allele_tuple = (locus, mut)
		# check for uniqueness
		if(allele_tuple in observed_alleles):
			continue
		cat['all'] += 1
		observed_alleles.append(allele_tuple)
		m = re.search(r'mbc\t(\w)\t(\d)\t(\w)\t(\d)\t(\w)\t(\d)\t(\w)\t(\d)', li)
		mbc_dict = dict()
		if m is not None:
			mbc_dict[m.group(1)] = int(m.group(2))
			mbc_dict[m.group(3)] = int(m.group(4))
			mbc_dict[m.group(5)] = int(m.group(6))
			mbc_dict[m.group(7)] = int(m.group(8))
		m = re.search(r'pbc\t(\w)\t(\d)\t(\w)\t(\d)\t(\w)\t(\d)\t(\w)\t(\d)', li)
		pbc_dict = dict()
		if m is not None:
			pbc_dict[m.group(1)] = int(m.group(2))
			pbc_dict[m.group(3)] = int(m.group(4))
			pbc_dict[m.group(5)] = int(m.group(6))
			pbc_dict[m.group(7)] = int(m.group(8))
		mdepth = 0
		pdepth = 0
		for base in mbc_dict:
			mdepth += mbc_dict[base]
		for base in pbc_dict:
			pdepth += pbc_dict[base]
		mbcfrac = dict()
		pbcfrac = dict()
		for base in mbc_dict:
			mbcfrac = float(mbc_dict[base])/float(mdepth)
		for base in pbc_dict:
			pbcfrac = float(pbc_dict[base])/float(pdepth)
		# is this actual discordance?
		# are we homoplasmic for both?
		if(mdepth > 9):
			cat['all_above_depth'].add(allele_tuple)
			mline = annotate_allele(allele_tuple,
					gb_frequencies, mitobank, phylo, veps, cat)
		if(pdepth > 9):
			cat['n_all_above_depth'].add(allele_tuple)
			pline = annotate_allele(allele_tuple,
					gb_frequencies, mitobank, phylo, veps, cat)
		max_shared_fraction = 0
		for base in mbc_dict:
			shared = min(pbcfrac, mbcfrac)
			if(shared > max_shared_fraction):
				max_shared_fraction = shared
		if((mdepth > 9) and (pdepth > 9)):
			if(max_shared_fraction <= 0.2):
				cat['almost_homoplasmy_change'].add(allele_tuple)
			if(max_shared_fraction == 0):
				cat['homoplasmy_change'].add(allele_tuple)


def annotate_allele(allele_tuple,
	gb_frequencies, mitobank, phylo, veps, cat,
	pathogenic_counts=dict(),
	fraction=0, individual='anon'):
	allele = allele_tuple
	not_in_mitobank = False
	cat['all'].add(allele)
	if(allele in phylo):
		cat['phylotree'].add(allele) 
	else:
		if(allele in mitobank):
			cat['mitobank'].add(allele)
			if allele in pathogenic_counts:
				pathogenic_counts[allele].append((individual, fraction))
			else:
				pathogenic_counts[allele] = [(individual, fraction)]
			if(is_rare(allele, gb_frequencies)):
				cat['mitobank_rare'].add(allele)
			else:
				cat['mitobank_common'].add(allele)
		else:
			not_in_mitobank = True
	if(not_in_mitobank):
		if(allele in gb_frequencies):
			observ = gb_frequencies[allele]
		else:
			observ = 0
		if(observ > 0):
			cat['not_mitobank_but_observed'].add(allele)
			if(allele in veps):
				siftline = veps[allele]
				siftdel = re.search(r'SIFT=deleterious', siftline)
				if(siftdel is not None):
					cat['not_mitobank_but_observed_deleterious'].add(allele)
				else:
					cat['not_mitobank_but_observed_nondeleterious'].add(allele)
			else:
				cat['not_mitobank_but_observed_nondeleterious'].add(allele)
			if(is_rare(allele, gb_frequencies)):
				cat['not_mitobank_but_observed_and_rare'].add(allele)
			else:
				cat['not_mitobank_but_observed_and_common'].add(allele)
		else:
			cat['novel'].add(allele)
			if(allele in veps):
				siftline = veps[allele]
				siftdel = re.search(r'SIFT=deleterious', siftline)
				if(siftdel is not None):
					cat['novel_deleterious'].add(allele)
				else:
					cat['novel_nondeleterious'].add(allele)
			else:
				cat['novel_nondeleterious'].add(allele)


def read_variant_freq_format(filename,
	gb_frequencies, mitobank, phylo, veps, cat,
	observed_alleles=[], pathogenic_counts=dict()):
	f = open(filename, 'r')
	num = 0
	for li in f:
		num += 1
		if(re.match("\A(\s*)$", li)):
			continue
		s = li.rstrip('\n').split('\t')
		m = re.match(r'\A(\d+)([a-zA-Z]+)', s[0])
		if m is not None:
			locus = int(m.group(1))
			mut = m.group(2)
			allele_tuple = (locus, mut)
			cat['all'] += 1
			observed_alleles.append(allele_tuple)
			number_of_ind = int(s[1])
			for index in range(1, number_of_ind + 1):
				annotate_allele(allele_tuple,
					gb_frequencies, mitobank, phylo, veps, cat,
					pathogenic_counts=pathogenic_counts,
					fraction=1,
					individual=str(num))
		else:
			print "Unrecognized allele format in %s" % li
			continue
	f.close()


def reset_cat():
	cat = {'all': set(),
	'all_above_depth': set(),
	'phylotree': set(),
	'mitobank': set(),
	'mitobank_rare': set(),
	'mitobank_common': set(),
	'not_mitobank_but_observed': set(),
	'not_mitobank_but_observed_nondeleterious': set(),
	'not_mitobank_but_observed_deleterious': set(),
	'not_mitobank_but_observed_and_common': set(),
	'not_mitobank_but_observed_and_rare': set(),
	'novel': set(),
	'novel_nondeleterious': set(),
	'novel_deleterious': set(),
	'almost_homoplasmy_change': set(),
	'homoplasmy_change': set()}
	return cat


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
	variant_glob_str = in_folder + '/low/*.ann_variant_freq'
	concordance_checked_glob_str = in_folder + '/low/*.concordance_checked'
	hets_glob_str = in_folder + '/low/*.hets'
	gb_frequencies = hettools.getGenbankFrequencies()  # key by tuple
	mitobank = hettools.readMitobank()  # key by tuple
	phylo = hettools.read_phylotree()  # key by tuple
	cat = reset_cat()
	observed_alleles = list()
	pathogenic_counts = dict()
	print "\n******************\n"
	if len(glob.glob(hets_glob_str)) == 0:
		print "No files match %s." % hets_glob_str
		sys.exit()
	veps = process_vep_files(hets_glob_str)
	read_hets_folder(hets_glob_str,
			gb_frequencies, mitobank, phylo, veps, cat,
			observed_alleles, pathogenic_counts=pathogenic_counts)
	print_cat(cat, mitobank, veps)
	print_pathogenic_counts(cat, mitobank, veps, pathogenic_counts)
	print_obs_alleles(observed_alleles)
	cat = reset_cat()

# This next section concerns homoplasmy and is temporarily excluded.
# This is old, bad code.
old_code = """
	observed_alleles = list()
	pathogenic_counts = dict()
	print "\n******************\n"
	for filename in glob.glob(variant_glob_str):
		basename = os.path.basename(filename)
		if(not re.match(r'all', basename)):
			continue
		print "variant_freq filename:%s" % basename
		m = re.match("([a-zA-Z0-9]+)\.(.*)variant_freq$", basename)
		if m is not None:
			keyword = m.group(1)
			vepfileguess = "homoplasmy/%s.vcf.fullvep" % keyword
			print "guessing the vep file is %s " % vepfileguess
			if(os.path.exists(vepfileguess)):
				print "\tfound vep file"
				veps = hettools.readInterestingVepLines(vepfileguess)
			else:
				print "\tNo such vep file"
				veps = {}
		else:
			print "\tno good guess for vep file"
			veps = {}
		read_variant_freq_format(filename,
			gb_frequencies, mitobank, phylo, veps, cat,
			observed_alleles, pathogenic_counts=pathogenic_counts)
	print_cat(cat, mitobank, veps)
	print_pathogenic_counts(cat, mitobank, veps, pathogenic_counts, alt_format=True, write_ind=False)
	cat = reset_cat()
	observed_alleles = list()
	print "\n******************\n"
	for filename in glob.glob(concordance_checked_glob_str):
		basename = os.path.basename(filename)
		if(not re.match(r'all', basename)):
			continue
		print "Opening file %s..." % basename
		# can we find the right vep file?
		m = re.match("([a-zA-Z0-9]+)\.(.*)concordance_checked$", basename)
		if m is not None:
			keyword = m.group(1)
			vepfileguess = "%s_interesting.vep" % keyword
			if(os.path.exists(vepfileguess)):
				print "\tfound vep file"
				veps = hettools.readInterestingVepLines(vepfileguess)
			else:
				print "\tNo such vep file"
				veps = {}
		else:
			print "\tno good guess for vep file"
			veps = {}
		read_concordance_checked_format(filename,
			gb_frequencies, mitobank, phylo, veps, cat, observed_alleles)
	print_cat(cat, mitobank, veps)
"""
