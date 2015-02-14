"""
Analyzes a hets/ folder.

Takes a keyword and .path file as input.

Expects a hets/ folder in the working directory (not the source directory).

Outputs a variety (too many) of analytical files to ./hets/burden/.

Kind of a mess.

Example:
python het_burden/het_burden.py -k 1kgenomes -p het_burden/lib/1kgenomes.path

"""
import os
import hettools
import het_concordance_tools
import categorize_alleles
from operator import itemgetter
import sys

def burden_by_family_position(cfFams,
	keyword, cutoff,
	mitobank, phylo):
	n_in_probands = dict()
	n_in_mothers = dict()
	n_in_fathers = dict()
	het_probands = set()
	het_mothers = set()
	het_fathers = set()
	non_het_probands = set()
	non_het_mothers = set()
	non_het_fathers = set()
	patho_in_probands = dict()
	patho_in_mothers = dict()
	patho_in_fathers = dict()
	phylo_in_probands = dict()
	phylo_in_mothers = dict()
	phylo_in_fathers = dict()
	for k in cfFams:
		fam = cfFams[k]
		n_in_mothers[k] = []
		n_in_fathers[k] = []
		n_in_probands[k] = []
		patho_in_probands[k] = []
		patho_in_mothers[k] = []
		patho_in_fathers[k] = []
		phylo_in_probands[k] = []
		phylo_in_mothers[k] = []
		phylo_in_fathers[k] = []
		for locus in fam.loci:
			# The 2 read cutoff doesn't correspond to any of the
			# main cutoff levels (low -> 1 read each strand),
			# but passing the low cutoff ensures we pass this
			# requirement. 
			if('proband' in fam.loci[locus]):
				het_probands.add(cfFams[k].proband)
				proband_obs = fam.loci[locus]['proband'][1]
				pbc = sorted(proband_obs.items(),
					key=itemgetter(1))
				pdepth = 0
				for i in pbc:
					pdepth += float(i[1])
				if(pbc[-2][1] > 1):
					n_in_probands[k].append(float(pbc[-2][1])/pdepth)
					if((locus, pbc[-2][0]) in mitobank):
						patho_in_probands[k].append(float(pbc[-2][1])/pdepth)
					if((locus, pbc[-2][0]) in phylo):
						phylo_in_probands[k].append(float(pbc[-2][1])/pdepth)
			if('mother' in fam.loci[locus]):
				het_mothers.add(cfFams[k].mother)
				mother_obs = fam.loci[locus]['mother'][1]
				mbc = sorted(mother_obs.items(),
					key=itemgetter(1))
				mdepth = 0
				for i in mbc:
					mdepth += float(i[1])
				if(mbc[-2][1] > 1):
					n_in_mothers[k].append(float(mbc[-2][1])/mdepth)
					# is there a pathogenic allele here?
					if((locus, mbc[-2][0]) in mitobank):
						patho_in_mothers[k].append(float(mbc[-2][1])/mdepth)
					if((locus, mbc[-2][0]) in phylo):
						phylo_in_mothers[k].append(float(mbc[-2][1])/mdepth)
			# prob_obs is a dict of integer values for 
			# number of times a base is observed.
			if('father' in fam.loci[locus]):
				het_fathers.add(cfFams[k].father)
				father_obs = fam.loci[locus]['father'][1]
				fbc = sorted(father_obs.items(),
					key=itemgetter(1))
				fdepth = 0
				for i in fbc:
					fdepth += float(i[1])
				if(fbc[-2][1] > 1):
					n_in_fathers[k].append(float(fbc[-2][1])/fdepth)
					if((locus, fbc[-2][0]) in mitobank):
						patho_in_fathers[k].append(float(fbc[-2][1])/fdepth)
					if((locus, fbc[-2][0]) in phylo):
						phylo_in_fathers[k].append(float(fbc[-2][1])/fdepth)
	all_probands = set()
	all_mothers = set()
	all_fathers = set()
	for k in cfFams:
		all_probands.add(cfFams[k].proband)
		all_mothers.add(cfFams[k].mother)
		all_fathers.add(cfFams[k].father)
		if(cfFams[k].proband not in het_probands):
			non_het_probands.add(cfFams[k].proband)
		if(cfFams[k].mother not in het_mothers):
			non_het_mothers.add(cfFams[k].mother)
		if(cfFams[k].father not in het_fathers):
			non_het_fathers.add(cfFams[k].father)
	if(not (os.path.exists("./hets/burden/"))):
		os.system("mkdir ./hets/burden")
	if(not (os.path.exists("./hets/burden/%s" % cutoff))):
		os.system("mkdir hets/burden/%s" % cutoff)
	total_p = len(list(het_probands))+ len(list(non_het_probands))
	total_m = len(list(het_mothers)) + len(list(non_het_mothers))
	total_f = len(list(het_fathers)) + len(list(non_het_fathers))
	base = "hets/burden/%s/%s.%s" % (cutoff, keyword, cutoff)
	li = """Individuals in cfFams:
		Probands: %i
		Mothers: %i
		Fathers: %i
		""" % (len(list(all_probands)),
			len(list(all_mothers)),
			len(list(all_fathers)))
	li += """Het burden by family position:
		het_probands/non_het: %i/%i\tsum=%i
		het_mothers/non_het: %i/%i\tsum=%i
		het_fathers/non_het: %i/%i\tsum=%i
		total individuals: %i
		""" % (
			len(list(het_probands)), len(list(non_het_probands)),
			total_p,
			len(list(het_mothers)), len(list(non_het_mothers)),
			total_m,
			len(list(het_fathers)), len(list(non_het_fathers)),
			total_f,
			total_p + total_m + total_f,
			)
	print li
	with open("%s.tmp_combined" % base, 'w') as outf:
		num_fam = 0
		li = "fam_num"
			# \tproband\tmother\tfather"
		li += "\tproband_depth\tmother_depth\tfather_depth"
		li += "\tproband_num_hets\tmother_num_hets\tfather_num_hets"
		li += "\tproband_max_het\tmother_max_het\tfather_max_het"
#		li += "\tproband_frac_phlyo\tmother_frac_phlyo\tfather_frac_phylo"
#		li += "\tproband_frac_patho\tmother_frac_patho\tfather_frac_patho"
		li += "\n"
		outf.write(li)
		for fam in cfFams:
			num_fam += 1
			afam = cfFams[fam]
			li = "%i" % num_fam
#			li += "\t%s\t%s\t%s" % (
#				afam.proband,
#				afam.mother,
#				afam.father)
			try:
				li += "\t%f\t%f\t%f" % (
					afam.depth['proband'],
					afam.depth['mother'],
					afam.depth['father'])
			except: 
				li += "\t0\t0\t0"
			li += "\t%f\t%f\t%f" % (
				len(n_in_probands[fam]),
				len(n_in_mothers[fam]),
				len(n_in_fathers[fam]))
			if(len(n_in_probands[fam]) > 0):
				highest_p = max(n_in_probands[fam])
			else:
				highest_p = 0
			if(len(n_in_mothers[fam]) > 0):
				highest_m = max(n_in_mothers[fam])
			else:
				highest_m = 0
			if(len(n_in_fathers[fam]) > 0):
				highest_f = max(n_in_fathers[fam])
			else:
				highest_f = 0
			li += "\t%f\t%f\t%f\n" % (
				highest_p,
				highest_m,
				highest_f)
			outf.write(li)
	with open("%s.proband_het_burden_by_num" % base, 'w') as outf:
		for fam in n_in_probands:
			num_het = len(n_in_probands[fam])
			outf.write("%i\n" % num_het)
	with open("%s.mother_het_burden_by_num" % base, 'w') as outf:
		for fam in n_in_mothers:
			num_het = len(n_in_mothers[fam])
			outf.write("%i\n" % num_het)
	with open("%s.father_het_burden_by_num" % base, 'w') as outf:
		for fam in n_in_fathers:
			num_het = len(n_in_fathers[fam])
			outf.write("%i\n" % num_het)
	with open("%s.proband_het_burden_by_highest_freq" % base, 'w') as outf:
		for fam in n_in_probands:
			if(len(n_in_probands[fam]) > 0):
				highest = max(n_in_probands[fam])
			else:
				highest = 0
			outf.write("%f\n" % highest)
	with open("%s.mother_het_burden_by_highest_freq" % base, 'w') as outf:
		for fam in n_in_mothers:
			if(len(n_in_mothers[fam]) > 0):
				highest = max(n_in_mothers[fam])
			else:
				highest = 0
			outf.write("%f\n" % highest)
	with open("%s.father_het_burden_by_highest_freq" % base, 'w') as outf:
		for fam in n_in_fathers:
			if(len(n_in_fathers[fam]) > 0):
				highest = max(n_in_fathers[fam])
			else:
				highest = 0
			outf.write("%f\n" % highest)
	with open("%s.proband_het_burden_by_pathogen" % base, 'w') as outf:
		for fam in patho_in_probands:
			num_het = len(patho_in_probands[fam])
			outf.write("%i\n" % num_het)
	with open("%s.mother_het_burden_by_pathogen" % base, 'w') as outf:
		for fam in patho_in_mothers:
			num_het = len(patho_in_mothers[fam])
			outf.write("%i\n" % num_het)
	with open("%s.father_het_burden_by_pathogen" % base, 'w') as outf:
		for fam in patho_in_fathers:
			num_het = len(patho_in_fathers[fam])
			outf.write("%i\n" % num_het)
	with open("%s.proband_het_burden_by_phylo" % base, 'w') as outf:
		for fam in phylo_in_probands:
			num_het = len(phylo_in_probands[fam])
			outf.write("%i\n" % num_het)
	with open("%s.mother_het_burden_by_phylo" % base, 'w') as outf:
		for fam in phylo_in_mothers:
			num_het = len(phylo_in_mothers[fam])
			outf.write("%i\n" % num_het)
	with open("%s.father_het_burden_by_phylo" % base, 'w') as outf:
		for fam in phylo_in_fathers:
			num_het = len(phylo_in_fathers[fam])
			outf.write("%i\n" % num_het)

def check_stats(cfFams):
	het_f = set()
	het_m = set()
	het_p = set()
	het_f_pos = set()
	het_m_pos = set()
	het_p_pos = set()
	for k in cfFams:
		for locus in cfFams[k].loci:
			alocus = cfFams[k].loci[locus]
			if('proband' in alocus):
				pbc = sorted(alocus['proband'][1].items(),
					key=itemgetter(1))
				het_p.add(cfFams[k].proband)
				if(pbc[-2][1] > 0):
					het_p_pos.add(cfFams[k].proband)
			if('mother' in alocus):
				het_m.add(cfFams[k].mother)
				mbc = sorted(alocus['mother'][1].items(),
					key=itemgetter(1))
				if(mbc[-2][1] > 0):
					het_m_pos.add(cfFams[k].mother)
			if('father' in alocus):
				het_f.add(cfFams[k].father)
				fbc = sorted(alocus['father'][1].items(),
					key=itemgetter(1))
				if(fbc[-2][1] > 0):
					het_f_pos.add(cfFams[k].father)
	li = """check_stats():
	het probands/positive freq hets: %i/%i
	het mothers/positive freq hets: %i/%i
	het fathers/positive freq hets: %i/%i
	""" % (
	len(list(het_p)), len(list(het_p_pos)),
	len(list(het_m)), len(list(het_m_pos)),
	len(list(het_f)), len(list(het_f_pos)),
	)
	print li


if __name__ == '__main__':
	src_path = os.path.dirname(os.path.realpath(__file__))
	args = het_concordance_tools.parse_input()
	keyword = 'all'
	# Relate individuals to families and files.
	paths_list = hettools.load_paths(args.paths, 'all', getAll=True)
	# Need the rcrs.
	refseq = hettools.read_crs(args.lib + '/rcrs.fa')
	# The following three values are universal library files,
	# all dicts with key=tuple of allele, used for annotation.
	gb_frequencies = hettools.getGenbankFrequencies()
	mitobank = hettools.readMitobank()
	phylo = hettools.read_phylotree()
	combined_cfFams = {
		'high': {},
		'mid': {},
		'low': {}}
	add_depths = False
	for path in paths_list:
		(idToFilename, bcFilenameToId) = hettools.read_map_file(path['map_file'])
		(probands, mothers, fathers, families, unused_ids
			) = hettools.read_in_stableid_families(path['trios_file'])
		# Read in the het file and fill in some cfHet objects.
		print "path=%s" % path['keyword']
		for cutoff in ['high', 'mid', 'low']:
			print "cutoff=%s" % cutoff
			hetsfile = "hets/%s/%s.%s.hets" % (cutoff, path['keyword'], cutoff)
			if(cutoff == 'high'):
				hetsfile = "hets/%s/%s.hets" % (cutoff, path['keyword'])
			hets = het_concordance_tools.read_hets_file(hetsfile, bcFilenameToId)

			cfFams = het_concordance_tools.convert_to_cfHet(hets, 
				probands, mothers, fathers, families)
			cfFams = het_concordance_tools.add_non_het_families(cfFams, families)
			#het_concordance_tools.analyze_cfFams(path, cfFams, idToFilename)
			if(add_depths):
				for fam in cfFams:
					cfFams[fam].add_depths(path['bc_folder'], idToFilename)
			combined_cfFams[cutoff].update(cfFams)
			print "Current cfFams stats:"
			check_stats(cfFams)
			print "Current combined cfFams stats:"
			check_stats(combined_cfFams[cutoff])
		#het_concordance_tools.write_concordance_files(path['hets_file'], cfFams, 
		#	gb_frequencies, refseq)

	# it is now possible to compare heteroplasmy burden
	#het_concordance_tools.write_concordance_files('hets/all', combined_cfFams, 
	#		gb_frequencies, refseq)
	for cutoff in ['high', 'mid', 'low']:
		het_loci = set()
		pm_pair = 0
		not_pm_pair = 0
		pm_above_cutoff = 0
		for k in combined_cfFams[cutoff]:
			[some_het_loci, _pm_pair, _not_pm_pair, _pm_above_cutoff
				] = combined_cfFams[cutoff][k].count_pairs()
			het_loci = het_loci.union(some_het_loci)
			pm_pair += _pm_pair
			not_pm_pair += _not_pm_pair
			pm_above_cutoff += _pm_above_cutoff
		li = """At cutoff %s:
		Families %i
		Heteroplasmic loci %i
		Proband-mother both heteroplasmic at loci %i
		Proband-mother both heteroplasmic at loci, with both above depth 10 %i
		Heteroplasmies not heteroplasmic in both proband and mother %i
		""" % (
		cutoff,
		len(combined_cfFams[cutoff]),
		len(list(het_loci)),
		pm_pair,
		pm_above_cutoff,
		not_pm_pair)
		print li
		burden_by_family_position(combined_cfFams[cutoff], 'all', cutoff, 
						mitobank, phylo)
