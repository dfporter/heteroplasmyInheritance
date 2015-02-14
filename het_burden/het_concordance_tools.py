"""
Contains a variety of utility functions for analyzing .hets file.

If run as __main__, creates .concordance.all and
.concordance.delterious files in 
hets/burden/cutoff/ folders from an input hets/ folder.

This program is an awful mess.

Example:

python het_burden/het_concordance_tools.py -p het_burden/lib/1kgenomes.path -k 1kgenomes

"""
import glob
import subprocess
import os
import re
from operator import itemgetter
import sys
import math
import argparse
import hettools
from cfHet import cfHet


def parse_input():
	src_path = os.path.dirname(os.path.realpath(__file__))
	parser = argparse.ArgumentParser(description=
	"""Process files indicating positions of heteroplasmy to generate concordance 
	and pathogenicity information.
	""")
	parser.add_argument('-k', '--keyword', default='none',
	help='(Required) Keyword')
	parser.add_argument('-l', '--lib', default="%s/lib/" % src_path,
	help="""Library folder contianing map file, trios file, rcrs.fa, 
	mitobank_variant_freqs.txt, and whatever other information files 
	are required. The filenames given for these files will be expected
	in the lib folder.
	""")
	parser.add_argument('--justVcf', dest='justVcf', 
	action='store_true', default=False, 
	help='Only parse .vcf files.')
	parser.add_argument('--noVep', dest='noVep', 
	action='store_true', default=False, 
	help='Assume hetsWithVep file has already been created and do not call vep.')
	parser.add_argument('--ignore_vep', dest='ignore_vep', 
	action='store_true', default=True, 
	help='Do not call vep and do not include its output.')
#	parser.add_argument('--by_cutoff', dest='by_cutoff', 
#	action='store_true', default=True, 
#	help='Write .concordance files by cutoff.')
	parser.add_argument('--use_filled', dest='use_filled', 
	action='store_true', default=False, 
	help='Use existing filled_in.hets files')
	parser.add_argument('-p', '--paths', default="%s/lib/paths2.path" % src_path,
	help="File of trios, map, bam folder, bc folder, ect. locations.")
	args = parser.parse_args() #args.justVcf now holds value
	print "Command line options given were:"
	return args


def read_hets_file(filename, bcFilenameToId, return_as_bc=False):
	"""For a given hetsbyindividual file, populate the global hets dict.
	filename -- hets file
	"""
	hets = dict()
	hetsF = open(filename, 'r')
	num_het_loci = 0
	missing_bc_files = set()
	found_bc_files = set()
	bc_files_found_by_guessing = set()
	for li in hetsF:
		num_het_loci += 1
		s = li.rstrip('\n').split('\t')
		bcFile =  s[0] + ".bc"
		locus = int(s[1])
		# s[0] contains no directory
		if(bcFile in bcFilenameToId):
			stable_id = bcFilenameToId[bcFile]
			hets.setdefault(stable_id, {locus: li})[locus] = li
			found_bc_files.add(bcFile)
		else:
			# Try guessing.
			trial_filename = "./" + bcFile
			if(os.path.exists(trial_filename)):
				m = re.match('(w+)\.mt\.bam\.bc', bcFile)
				if(m is not None):
					stable_id = m.group(1)
					hets.setdefault(stable_id, {locus: li})[locus] = li
					bcFilenameToId[bcFile] = stable_id
					bc_files_found_by_guessing.add(bcFile)
					found_bc_files.add(bcFile)
				else:
					missing_bc_files.add(bcFile)
			else:
				missing_bc_files.add(bcFile)
	#	print "adding %s id to hets dict from -%s- filename in hets file" % (stable_id, bamFile)
	hetsF.close()
	num_stable_ids_with_het_loci = len(found_bc_files) + len(missing_bc_files)
	print """read_hets_file()
		het loci %i
		number of stable ids with het loci %i
		bc files found %i
		bc files found only by guessing %i
		bc files we could not find %i
		The bc files we could not find are: %s\n""" % (
		num_het_loci, 
		len(found_bc_files) + len(missing_bc_files),
		len(found_bc_files), 
		len(bc_files_found_by_guessing),
		len(missing_bc_files),
		str(missing_bc_files))
	if(return_as_bc):
		for stable_id in hets:
			for locus in hets[stable_id]:
				bcre = re.search(
r'\[\(\'(\w)\', (\d+)\), \(\'(\w)\', (\d+)\), \(\'(\w)\', (\d+)\), \(\'(\w)\', (\d+)\)', 
					hets[stable_id][locus]).groups()
				bc = {bcre[0]: int(bcre[1]), bcre[2]: int(bcre[3]),
					bcre[4]: int(bcre[5]), bcre[6]: int(bcre[7])}
				hets[stable_id][locus] = sorted(bc.items(), key=itemgetter(1))
	return hets


def convert_to_cfHet(hets, probands, mothers, fathers, families):
	"""Convert the hets dict to a dict of cfHets objects.
	The cfHets objects are organized into a cfFams dict with 
	key=(proband, mother) tuple, value=cfHet object.
	The various het loci are added to the respective cfHet objects
	using the class method .set_locus_from_hets_file.
	The purpose of the dict of cfHet objects is to be able to
	write concordance files, using the cfHet methods to deal with
	organizational issues such as allele output order.
	"""
	cfFams = dict()  # key= (proband,mother) tuple. value=a cfHet object. 
	num_have_vep_annotation = 0
#	print "we are aware of the sets: probands:%s,\n\nmothers:%s\n\nfathers:%s\n\n" % (
#		str(probands), str(mothers), str(fathers) )
	het_probands = set()
	het_mothers = set()
	het_fathers = set()	
	did_not_find_family = set()
	for stable_id in hets:
		for locus in hets[stable_id]:
			vep = ""
			s = hets[stable_id][locus].rstrip('\n').split('\t')
			maf = s[2]
			vepSearch = re.search(r'\t\$([^\$]+)', hets[stable_id][locus])
			if (vepSearch is not None):
				vep = vepSearch.group(1)
			bcre = re.search(r'\[\(\'(\w)\', (\d+)\), \(\'(\w)\', (\d+)\), \(\'(\w)\', (\d+)\), \(\'(\w)\', (\d+)\)', hets[stable_id][locus] ).groups()
			bc = {bcre[0]: bcre[1], bcre[2]: bcre[3],
				bcre[4]: bcre[5], bcre[6]: bcre[7]}
			found_in_a_family = False
			if(stable_id in probands):
				found_in_a_family = True
				het_probands.add(stable_id)
				proband = stable_id
				mother = families[probands[stable_id]][1]
				father = families[probands[stable_id]][2]
				add_from_hets_file(cfFams, proband, mother, father, 
							stable_id, locus, maf, bc, vep)
			if(stable_id in mothers):
				found_in_a_family = True
				het_mothers.add(stable_id)
				mother = stable_id
				proband = families[mothers[stable_id]][0]
				father = families[mothers[stable_id]][2]
				add_from_hets_file(cfFams, proband, mother, father, 
							stable_id, locus, maf, bc, vep)
			if(stable_id in fathers):
				found_in_a_family = True
				het_fathers.add(stable_id)
				father = stable_id
				proband = families[fathers[stable_id]][0]
				mother = families[fathers[stable_id]][1]
				add_from_hets_file(cfFams, proband, mother, father, 
							stable_id, locus, maf, bc, vep)
			if(not found_in_a_family):
				did_not_find_family.add(stable_id)
	print """convert_to_cfHet():
		Individuals: %i 
		Heteroplasmic probands: %i 
		Heteroplasmic mothers: %i 
		Heteroplasmic fathers: %i 
		Individuals we could not find in a family: %i\n""" % (
		len(hets),
		len(het_probands),
		len(het_mothers),
		len(het_fathers),
		len(did_not_find_family))	
	return cfFams


def add_from_hets_file(cfFams, proband, mother, father, stable_id, locus, maf, bc, vep):
	k = (proband, mother)
	if(k in cfFams):
		cfFams[k].set_locus_from_hets_file(
			locus, stable_id, maf, bc, vep)
	else:  # new family
		cfFams[k] = cfHet(proband=proband, mother=mother, father=father)
		cfFams[k].set_locus_from_hets_file(
			locus, stable_id, maf, bc, vep)


def analyze_cfFams(path, cfFams, idToFilename, cutoff='high'):
	"""Fills in the info in cfHet objects.
	The path argument is a path dict, used to retrieve a bc_folder.
	The cfFams dict passed as argument holds key=(proband, mother) and
	value=cfHet object. The cfHet objects have loci set from the hets file,
	and therefore miss situations where one individual is not in the hets
	file. We need to be able to compare all mother/proband pairs,
	so this function adds info from the .bc files.
	Currently, this function does a depth check on the added individual.
	It requires the pair to be at least 10 total reads deep. This is a less
	stringent cutoff than the lowest heteroplasmy cutoff (minimum 5 reads 
	each direction), so it will never exclude something that would be
	included by the findHets program.
	"""
	both_individuals_in_hets_file = 0
	only_one_in_hets_file = 0
	missing_individuals = set()
	properties_of_pair = {
	'pair_above_10_per_strand_depth': 0,
	'pair_above_10_total_depth': 0,
	'pair_below_total_depth': 0,
	'missing_locus': set()}
	for k in cfFams:
		for locus in cfFams[k].loci:
			#are we paired?
			if(('mother' in cfFams[k].loci[locus])
				and ('proband' in cfFams[k].loci[locus])):
				# Paired
				both_individuals_in_hets_file += 1
			else:
				only_one_in_hets_file += 1
	for k in cfFams:
		for locus in cfFams[k].loci:
			# avoid double-counting at a pair
			dont_count_properties=False
			#look up the .bc file
			bcfolder = path['bc_folder'].rstrip('/') + '/'
			if('proband' not in cfFams[k].loci[locus]):
				unknown_individual = cfFams[k].proband
				if(unknown_individual not in idToFilename):
					missing_individuals.add(unknown_individual)
					cfFams[k].set_locus_from_bc_file(locus, unknown_individual, 0, 
					'Error: could not find file', 
					{'A': 0, 'T': 0, 'C': 0, 'G': 0})
				else:
					bcfilename = bcfolder + idToFilename[unknown_individual]
					add_locus_from_bc(cfFams[k], unknown_individual,
						bcfilename, locus, properties_of_pair,
						dont_count_properties=dont_count_properties,
						cutoff=cutoff)
					dont_count_properties=True
			if('mother' not in cfFams[k].loci[locus]):
				unknown_individual = cfFams[k].mother
				if(unknown_individual not in idToFilename):
					missing_individuals.add(unknown_individual)
					cfFams[k].set_locus_from_bc_file(locus, unknown_individual, 0, 
					'Error: could not find file', 
					{'A': 0, 'T': 0, 'C': 0, 'G': 0})
				else:
					bcfilename = bcfolder + idToFilename[unknown_individual]
					add_locus_from_bc(cfFams[k], unknown_individual,
						bcfilename, locus, properties_of_pair,
						dont_count_properties=dont_count_properties,
						cutoff=cutoff)
					dont_count_properties=True
			if('father' not in cfFams[k].loci[locus]):
				unknown_individual = cfFams[k].father
				if(unknown_individual not in idToFilename):
					print "Could not find %s" % unknown_individual
					print "Looked in file %s" % bcfilename
					missing_individuals.add(unknown_individual)
					cfFams[k].set_locus_from_bc_file(locus, unknown_individual, 0, 
					'Error: could not find file', 
					{'A': 0, 'T': 0, 'C': 0, 'G': 0})
				else:
					bcfilename = bcfolder + idToFilename[unknown_individual]
					add_locus_from_bc(cfFams[k], unknown_individual,
						bcfilename, locus, properties_of_pair,
						dont_count_properties=True,
						cutoff=cutoff)
	print """analyze_cfFams(): 
		both individuals in hets file %i 
		only one in hets file %i 
		pair is above 10 total depth (loci can be compared) %i 
		missing individual %i 
		missing locus %i 
		pair is above 10 per strand depth %i 
		pair is below 10 total depth %i\n""" % (
		both_individuals_in_hets_file,
		only_one_in_hets_file,
		both_individuals_in_hets_file + properties_of_pair['pair_above_10_total_depth'],
		len(list(missing_individuals)),
		len(list(properties_of_pair['missing_locus'])),
		properties_of_pair['pair_above_10_per_strand_depth'],
		properties_of_pair['pair_below_total_depth'])


def add_locus_from_bc(cfHet_obj, unknown_individual, bcfilename,
			locus, properties_of_pair,
			dont_count_properties=False, cutoff='high'):
	"""cutoff is passed as an argument (holding a string like 'high')
	in case we use it later. For the moment, we just use 10 as a cutoff
	for everything
	"""
	bcfile = open(bcfilename, 'r')
	found_paired_locus = False
	cutoff = 10
	for li in bcfile: #chr locus ref_allele basecounts
		q = li.rstrip('\n').split('\t')
		if(int(q[1]) != locus):
			continue
		found_paired_locus = True
		# Is this locus usable and het?
		totaldepth = 0
		totalfor = 0
		totalrev = 0
		# Convert to int.
		s = q[0:3]
		for i in q[3:11]:
			s.append( int(i) )
		for i in s[3:7]:
			totalfor += i
		for i in s[7:11]:
			totalrev += i
		obsBases = dict()
		obsAlleles = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
		i = 3
		for a in ['A', 'T', 'C', 'G']:
			obsAlleles[a] = s[i] + s[4+i]
			i += 1
		if((totalfor + totalrev) >= cutoff):
			sorted_bc = sorted(obsAlleles.items(), key=itemgetter(1))
			maf = 1 - (float(sorted_bc[-1][1])/float(
				totalfor + totalrev))
			if(not dont_count_properties):
				properties_of_pair['pair_above_10_total_depth'] += 1	
			if((totalfor > 9) and (totalrev > 9)): 
				# This site is usable.
				cfHet_obj.set_locus_from_bc_file(
				locus, unknown_individual, maf, 
				'bc file above 10 reads each direction', obsAlleles)
				if(not dont_count_properties):
					properties_of_pair['pair_above_10_per_strand_depth'] += 1
			else:
				cfHet_obj.set_locus_from_bc_file(
				locus, unknown_individual, maf, 
				"bc file below 10for/10rev reads: %i depth" % (
				totalfor + totalrev), obsAlleles )
		else:
			cfHet_obj.set_locus_from_bc_file(
			locus, unknown_individual, 0, 
			'Error: low depth', obsAlleles)
			if(not dont_count_properties):
				properties_of_pair['pair_below_total_depth'] += 1
	if(not found_paired_locus):
		if(not dont_count_properties):
			properties_of_pair['missing_locus'].add(unknown_individual)
	bcfile.close()


def concordance_of_locus(k, locus, selection, 
				cfFams, gb_frequencies, refseq, 
				include_unused=False,
				mode='mitobank', mitobank={}):
	"""Return formatted lines of concordance information.
	"""
	if(selection=='all'):
		line = cfFams[k].write_frequencies(locus, 'all')
		majAllele = cfFams[k].alleles[locus][0]
		minAllele = cfFams[k].alleles[locus][1]
		if((locus, majAllele) in gb_frequencies):
			maj_freq = str(gb_frequencies[(locus, majAllele)])
		else:
			if(majAllele == refseq[locus-1]):
				maj_freq = '24000' # This is the reference allele.
			else:
				maj_freq = str(0)
		if((locus, minAllele) in gb_frequencies):
			min_freq = str(gb_frequencies[(locus, minAllele)])
		else:
			if(minAllele == refseq[locus-1]):
				min_freq = '24000' # This is the reference allele.
			else:
				min_freq = str(0)
		lineForPlot = str(locus) + "\t" 
		lineForPlot += "\t".join([majAllele, str(maj_freq), minAllele, str(min_freq), "\t"]) 
		lineForPlot += cfFams[k].write_frequencies(locus, 'all') 
		lineForPlot += "\t" + "\t".join([cfFams[k].mother, cfFams[k].proband])
		if(include_unused):
			family_number = probands[cfFams[k].proband]
			lineForPlot += "\t" + unused_ids[family_number] + "\n"
		else:
			lineForPlot += "\n"
		return (line, lineForPlot)
	if(selection=='deleterious'):  # deleterious
		line = cfFams[k].write_frequencies(locus, 'deleterious')
		majAllele = cfFams[k].alleles[locus][0]
		minAllele = cfFams[k].alleles[locus][1]
		if((locus, majAllele) in gb_frequencies):
			maj_freq = str(gb_frequencies[(locus, majAllele)])
		else:
			if(majAllele == refseq[locus-1]):
				maj_freq = '24000' # this is the reference allele
			else:
				maj_freq = str(0)
		if((locus, minAllele) in gb_frequencies):
			min_freq = str(gb_frequencies[(locus, minAllele)])
		else:
			if(minAllele == refseq[locus-1]):
				min_freq = '24000' # this is the reference allele
			else:
				min_freq = str(0)
		lineForPlot = str(locus) + "\t" 
		lineForPlot += "\t".join([majAllele, str(maj_freq), minAllele, str(min_freq), "\t"]) 
		lineForPlot += cfFams[k].write_frequencies(locus, 'deleterious',
				mode='mitobank', mitobank=mitobank) 
		lineForPlot += "\t" + "\t".join([cfFams[k].mother, cfFams[k].proband])
		# unused id
		if(include_unused):
			family_number = probands[ cfFams[k].proband ]
			lineForPlot += "\t" + unused_ids[family_number] + "\n"
		else:
			lineForPlot += "\n"
		return ("", lineForPlot)


def add_non_het_families(cfFams, families):
	"""Adds extra families from given dict to cfFams obj.
	Returns the updated cfFams object.
	"""
	known_ind = set()
	known_fams = set()
	ind_in_trio_not_in_hets = set()
	fam_in_trio_not_in_hets = set()
	for k in cfFams:
		known_ind.add(cfFams[k].proband)
		known_ind.add(cfFams[k].mother)
		known_ind.add(cfFams[k].father)
		known_fams.add((cfFams[k].proband, cfFams[k].mother))
	for atrio in families:
		proband = families[atrio][0]
		mother = families[atrio][1]
		father = families[atrio][2]
		for person in [x for x in families[atrio][0:3] if x not in known_ind]:
			ind_in_trio_not_in_hets.add(person)
		if((proband, mother) not in known_fams):
			k = (proband, mother)
			fam_in_trio_not_in_hets.add(k)
			cfFams[k] = cfHet(proband=proband, mother=mother, father=father)
	li = """add_non_het_families():
	individuals known from het families %i
	individuals in trios file, but not in hets file %i
		...total individuals %i
	families from trios file %i
	families known from hets file %i
	families unknown from hets file %i
		...total families now included in cfFams %i
	""" % (
	len(list(known_ind)),
	len(list(ind_in_trio_not_in_hets)),
	len(list(known_ind)) + len(list(ind_in_trio_not_in_hets)),
	len(families),
	len(list(known_fams)),
	len(list(fam_in_trio_not_in_hets)),
	len(cfFams))
	print li
	return cfFams


def write_concordance_files(file_basename, keyword, hetsfile, cfFams, 
				gb_frequencies, refseq):
	"""Writes .concordance.all and .concordance.deleterious files.
	Takes a filename (to set the output filenames) and a cfFams 
	object as arguments.
	"""
	basename_output = file_basename
	basename_output = re.sub(r'\.hetsByIndividual','', basename_output) 
	print "write_concordance_files():"
	print "\toutput basename for concordance files: %s" % basename_output
	freq_all = open("%s.concordance.all" % basename_output, 'w')
	freq_bad = open("%s.concordance.deleterious" % basename_output, 'w')
	mitobank = hettools.readMitobank()
	for k in cfFams:
		line = "" # For output to terminal.
		for locus in cfFams[k].concord:
			locus = int(locus)
			line += cfFams[k].write_concordance(locus)	
			ok_to_output = cfFams[k].set_frequencies(locus)
			if(not ok_to_output):
				continue
			line += "\t"
			if(cfFams[k].write_frequencies(locus, 'all') and 
				cfFams[k].is_above_depth(locus)): 
				# Any pathogenicity or lack thereof.
				(lineExt, lineForPlot) = concordance_of_locus(k, locus, 'all', 
					cfFams, gb_frequencies, refseq)
				line += lineExt
				freq_all.write(lineForPlot)
			if(cfFams[k].write_frequencies(locus, 'deleterious',
					mode='mitobank', mitobank=mitobank) and
				cfFams[k].is_above_depth(locus)): 
				# Deleterious.
				(lineExt, lineForPlot) = concordance_of_locus(k, locus, 
				'deleterious', cfFams, gb_frequencies, refseq,
				mode='mitobank', mitobank=mitobank)
				freq_bad.write(lineForPlot)
	freq_all.close()
	freq_bad.close()


def write_filled_in_hets(cfFams, outfile):
	"""Writes each family in a cfFams object to a file.
	"""
	outF = open(outfile, 'w')
	for fam in cfFams:
		li = cfFams[fam].write_hets_format()
		outF.write(li)
	outF.close()


def process_a_path(hetsfilename, path, keyword, args, refseq, cutoff):
	"""Creates a cfFams object from a hets filename.
	Returns the cfFams object.
	"""
	print "keyword=%s" % keyword
	print "hetsfile=%s" % hetsfilename
	(idToFilename, bcFilenameToId) = hettools.read_map_file(path['map_file'])
	(probands, mothers, fathers, families, unused_ids
		) = hettools.read_in_stableid_families(path['trios_file'])
	# Unless requested otherwise, call vep and add vep information
	# return value is the basename of path['hets_file'] + .vcf
	hets_with_vep_filename = vep_io_on_path(
				hetsfilename, keyword, args, refseq)
	hets_filename_to_read = hets_with_vep_filename
	if(args.ignore_vep):
		hets_filename_to_read = hetsfilename
	if(not args.justVcf):
		cfFams = build_cfFams_on_path(hets_filename_to_read,
			path, keyword, args, refseq,
			idToFilename, bcFilenameToId,
			probands, mothers, fathers, families,
			cutoff=cutoff)
	else:
		vep_io_on_path("hets/%s.hets" % keyword, keyword, args, refseq)
	return cfFams


def vep_io_on_path(hetsfile, keyword, args, refseq):
	"""Is this function used?"""
	vcf_filename = hettools.writeVcf(hetsfile, args.noVep, refseq) 
	if(not args.ignore_vep):
		hettools.update_with_vep_output(hetsfile, vcf_filename, 
			"%sWithVep.txt" % hetsfile , 
			keyword=keyword, 
			doCallVep=not(args.noVep))
	return "%sWithVep.txt" % hetsfile


def build_cfFams_on_path(hets_with_vep_filename, 
		path, keyword, args, refseq,
		idToFilename, bcFilenameToId,
		probands, mothers, fathers, families, 
		cutoff='high'):
	# Read lines of hets file into dict hets, with key=stable_id, value=file line.
	if(args.use_existing_filled_in_files):
		filled_in_fname = "hets/%s/%s.%s.filled_in.hets" % (cutoff, keyword, cutoff)
		hets = read_hets_file(filled_in_fname, bcFilenameToId)
		cfFams = convert_to_cfHet(hets, probands, mothers, fathers, families) 
		gb_ref_filename = args.lib.rstrip('/') + "/mitobank_variant_freqs.txt"
		gb_frequencies = hettools.getGenbankFrequencies(gb_ref_filename)
		write_concordance_files("hets/%s/%s.%s" % (cutoff, keyword, cutoff), 
				keyword, path['hets_file'], cfFams, 
				gb_frequencies, refseq)
		return cfFams
	hets = read_hets_file(hets_with_vep_filename, bcFilenameToId)
	# Read in all the data from the hets file (in hets object) into a cfFams dict, 
	# key=tuple of names, value=cfHet objects.
	cfFams = convert_to_cfHet(hets, probands, mothers, fathers, families) 
	analyze_cfFams(path, cfFams, idToFilename, cutoff=cutoff)
	# Genbank frequencies.
	gb_ref_filename = args.lib.rstrip('/') + "/mitobank_variant_freqs.txt"
	gb_frequencies = hettools.getGenbankFrequencies(gb_ref_filename)
	write_concordance_files("hets/%s/%s.%s" % (cutoff, keyword, cutoff), 
			keyword, path['hets_file'], cfFams, 
			gb_frequencies, refseq)
	write_filled_in_hets(cfFams, "hets/%s/%s.%s.filled_in.hets" % (cutoff, keyword, cutoff))
	return cfFams


def loop_through_files_at_cutoff(keyword, args, cutoff, do_loop=True):
	"""Top level function below __main__.
	For each keyword in the paths file, expects a hets/cutoff/.hets
	file to exists. Calls process_a_path on each .hets file and makes
	a combined cfFams object. Combines all the .concordance files
	by a cat command at the end.
	do_loop is not implemented yet. Assumed true.
	"""
	print "cutoff=%s" % cutoff
	if(not os.path.exists("hets/%s/" % cutoff)):
		print "Expect the folder hets/%s/ to exist... Exiting..." % cutoff
		sys.exit()
	paths_list = hettools.load_paths(args.paths, keyword, getAll=True)
	for path in paths_list:
		keyword=path['keyword']
		if(cutoff == 'high'):
			hetsfilename = "hets/%s/%s.hets" % (cutoff, keyword)
		else:
			hetsfilename = "hets/%s/%s.%s.hets" % (cutoff, keyword, cutoff)
		cfFams_of_path = process_a_path(hetsfilename, path, keyword, args, refseq,
					cutoff=cutoff)
		if(cfFams_of_path is not None):
			all_cfFams.update(cfFams_of_path)
	write_filled_in_hets(all_cfFams, 
		"hets/%s/%s.%s.filled_in.hets" % (cutoff, keyword, cutoff))
	cmdl_all = ["cat"]
	cmdl_bad = ["cat"]
	for path in paths_list:
		basename_ = "hets/%s/%s.%s" % (cutoff, path['keyword'], cutoff)
		basename_ = re.sub(r'\.hetsByIndividual','', basename_) 
		freq_all = "%s.concordance.all" % (basename_)
		freq_bad = "%s.concordance.deleterious" % (basename_)
		cmdl_all.append(freq_all)
		cmdl_bad.append(freq_bad)
	all_outfilename = "hets/%s/all.%s.concordance.all" % (cutoff, cutoff)
	bad_outfilename = "hets/%s/all.%s.concordance.deleterious" % (cutoff, cutoff)
	print "__main__: concatenating .concordance files"
	print "\t%s" % str(cmdl_all)
	with open(all_outfilename, 'w') as outfile:
		subP = subprocess.Popen(cmdl_all, stdout=outfile)
		subP.wait()
	with open(bad_outfilename, 'w') as outfile:
		subP = subprocess.Popen(cmdl_bad, stdout=outfile)
		subP.wait()


if __name__ == '__main__':
	args = parse_input()
	keyword = args.keyword
	refseq = hettools.read_crs(args.lib + '/rcrs.fa')
	all_cfFams = dict()
	args.use_existing_filled_in_files = args.use_filled
	print args.paths
	if(keyword=='all'):
#		if(args.by_cutoff):
		loop_through_files_at_cutoff(keyword, args, cutoff='high')
		loop_through_files_at_cutoff(keyword, args, cutoff='mid')
		loop_through_files_at_cutoff(keyword, args, cutoff='low')
#		else:
#			loop_through_files_at_cutoff(keyword, args, cutoff='high')
	else:
		path = hettools.load_paths(args.paths, keyword)
#		if(args.by_cutoff):
		loop_through_files_at_cutoff(keyword, args, cutoff='high',
			do_loop=False)
		loop_through_files_at_cutoff(keyword, args, cutoff='mid',
			do_loop=False)
		loop_through_files_at_cutoff(keyword, args, cutoff='low',
			do_loop=False)
#		else:
#			loop_through_files_at_cutoff(keyword, args, cutoff='high',
#				do_loop=False)
