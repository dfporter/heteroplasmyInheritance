import subprocess
import os
import re
import inspect
from operator import itemgetter


# =============================
# Utility functions for mtDNA analysis
# =============================


def load_paths(pathsfilename, keyword, getAll=False):
	src_path = os.path.dirname(os.path.realpath(__file__))
        foundKeyword = False
        f = open(pathsfilename, 'r')
	paths = dict()
	paths_list = []
	if(getAll):
		for li in f:
			if(re.match(r'\A#', li)):
				continue
			s = li.rstrip('\n').split('\t')
			foundKeyword = True
			paths_list.append({ 
				'keyword': s[0],
				'bam_folder': s[1],
				'bc_folder': s[2],
				'map_file': s[3],
				'trios_file': s[4],
				'hets_file': s[5],
				'vcf_file': s[6]})
	else:
		for li in f:
			if(re.match(r'\A#', li)):
				continue
			s = li.rstrip('\n').split('\t')
			if(s[0]==keyword):
				foundKeyword = True
				paths = {
					'keyword': s[0],
					'bam_folder': s[1],
					'bc_folder': s[2],
					'map_file': s[3],
					'trios_file': s[4],
					'hets_file': s[5],
					'vcf_file': s[6]}
        f.close()
	for path in paths_list:
		for col in path:
			path[col] = re.sub('src_path', src_path, path[col])
        if(foundKeyword):
		if(getAll):
			return paths_list
		else:
                	return paths
        else:
                return False


def is_bad_loci(locus):
        """Returns True if a locus is not usable, False otherwise.
        Bad loci are defined as those not passing a >10% cutoff.
        """
	badLoci = { 2487: 1, 3106: 1, 3107: 1, 3492: 1, 5734: 1, 8859: 1, 8860: 1,
            8861: 1, 8862: 1, 10306: 1 }
        locusI = int(locus)
        if( (locusI<71) or (locusI>16462) or ( (locusI<316) and (locusI>310) ) ):
                return True
        if(locusI in badLoci):
                return True
        return False


def catch_CRS(locus, genotype):
	"""Called to catch alleles that are just errors in the crs.
	"""
	rareInCRS =  'WT(CRS contains rare polymorphism)'
	errorInCRS = 'WT(CRS contains an error)'
	crs = {
	'A263G': rareInCRS,
	'TCCCCC310TCCCCCC': rareInCRS,
	'A750G': rareInCRS,
	'A1438G': rareInCRS,
	'CN3106C': errorInCRS,
	'G3423T': errorInCRS,
	'A4769G': rareInCRS,
	'G4985A': errorInCRS,
	'A8860G': rareInCRS,
	'G9559C': errorInCRS,
	'T11335C': errorInCRS,
	'G13702C': errorInCRS,
	'G14199T': errorInCRS,
	'G14272C': errorInCRS,
	'G14365C': errorInCRS,
	'G14368C': errorInCRS,
	'T14766C': errorInCRS,
	'A15326G': rareInCRS,
	}
	crs_lookup = dict()
	for mut in crs:
       	 	s = re.match(r'\A([AGCTN]+)([0-9]+)([AGCTN]+$)', mut)
	        crs_lookup[ (s.groups()[1], s.groups()[2]) ] = True
        if ((str(locus), genotype)  in crs_lookup):
                # This is an allele that is an error or rare polymorphism
		# in the CRS.
                return True
        else:
                return False


def read_map_file(mapFilename):
        """Reads in a file of <id \t filename>  (.map) information.
        The filename is the basename of the .bam.bc file conatining all reads.
	Two types of IDs are allowed in this project, person-IDs and
	file-IDs. Both IDs should have a unique 1:1 mapping with a 
	.bam/.bc file and an individual. There is therefore no functional
	difference.
	This funciton sets the variable:
       	idToFilename (dict) as key=file id, value=file basename.
	Crucial to these programs working is that the basecall filename
	is the bam file name with a .bc suffix, and in a different folder.
	The bam and bc folders are given as arguments, and the map filename
	is used to connect file basenames with file ids.
        """
        mapF = open(mapFilename, 'r')
	idToFilename = dict()
	bcFilenameToId = dict()
        for l in mapF:
                s = l.rstrip('\n').split('\t')
		basename = os.path.basename(s[1])
                idToFilename[s[0]] = basename  # File is <id \t file basename>
		bcFilenameToId[basename] = s[0]  # Key is by basename, not the full path
        mapF.close()
	return (idToFilename, bcFilenameToId) 


def read_all_map_files(paths_filename):
	paths_list = load_paths(paths_filename, 'all', getAll=True)
	full_bc_path_to_id = dict()
	id_to_full_filename = dict()
	for path in paths_list:
		(id_to_filename, bc_filename_to_id
			) = read_map_file(path['map_file'])
		for bc_filename in  bc_filename_to_id:
			full_path = path['bc_folder'] + '/' + bc_filename
			full_bc_path_to_id[full_path] = bc_filename_to_id[bc_filename]
		for stableid in id_to_filename:
			full_path = path['bc_folder'] + '/' + id_to_filename[stableid]
			id_to_full_filename[stableid] = full_path
	return (id_to_full_filename, full_bc_path_to_id)


def read_in_stableid_families(familyFilename):
        """Reads in a file of pedigree information.
        Sets the global variable families (dict) with file IDs.
	Stable IDs are person-IDs. These
	IDs are assumed to be redundant.
        File must be in the format:
        unused_id     proband_stable_id       mother_stable_id        father_stable_id
        gender  proband_path    proband_file_id       dad_path        dad_file_id   
        mother_path        mother_file_id
        Except only the file_id is used. 
        Proband id = s[6]. Father id = s[8]. Mother id = s[10]
        """
        trios_file = open(familyFilename,'r')
        num = 0
	probands = dict()
	mothers = dict()
	fathers = dict()
	families = dict()
	unused_ids = dict()
        for l in trios_file:
                num += 1
                s = l.rstrip('\n').split('\t')
		probands[s[1]] = num
		mothers[s[2]] = num
		fathers[s[3]] = num
		families[num] = [s[1], s[2], s[3]]
		unused_ids[num] = s[0]
        trios_file.close()
	return (probands, mothers, fathers, families, unused_ids)


def read_in_file_id_to_stable_id(trios_filename):
	file_to_stable = dict()
	trios_file = open(trios_filename,'r')
	file_trios = dict()
	trio_num = 0
	file_to_trio_num = dict()
	for li in trios_file:
                s = li.rstrip('\n').split('\t')
		file_to_stable[s[6]] = s[1]
		file_to_stable[s[10]] = s[2]
		file_to_stable[s[8]] = s[3]
		file_trios[trio_num] = [s[1], s[2], s[3]]
		file_to_trio_num[s[6]] = trio_num
		file_to_trio_num[s[8]] = trio_num
		file_to_trio_num[s[10]] = trio_num
		trio_num += 1
	trios_file.close()
	return file_to_stable


def read_in_families(familyFilename, idToFilename):
        """Reads in a file of pedigree information.
        Sets the global variable families (dict) with file IDs.
        File must be in the format:
        unused_id     proband_stable_id       mother_stable_id        father_stable_id
        gender  proband_path    proband_file_id       dad_path        dad_file_id   
        mum_path        mum_file_id
        Except only the file_id is used. 
        Proband id = s[6]. Father id = s[8]. Mother id = s[10]
        """
        trios_file = open(familyFilename,'r')
        num = 0
	probands = dict()
	mothers = dict()
	fathers = dict()
	families = dict()
	unused_ids = dict()
        for li in trios_file:
                num += 1
                s = li.rstrip('\n').split('\t')
                if(s[8] in idToFilename): 
			# idToFilename has value=full path.
			basename = os.path.basename(idToFilename[s[8]])
			probands[s[6]] = num
			mothers[s[10]] = num
			fathers[s[8]] = num
			families[num] = [s[6], s[10], s[8]]
			unused_ids[num] = s[0]
        trios_file.close()
	return (probands, mothers, fathers, families, unused_ids)


def read_crs(filename=r'lib/rcrs.fa'):
        """Read in the revised cambridge reference sequence.
        filename -- full path to rcrs.fasta
        Sets global variable refseq to be the rCRS.
        """
	src_path = os.path.dirname(os.path.realpath(__file__))
	filename = "%s/lib/rcrs.fa" % src_path
	refseq = ""
        rCRSf = open(filename, 'r')
        for li in rCRSf:
                if (re.match('>', li)):
                        continue
                li = li.rstrip('\n')
                refseq += li
        rCRSf.close()
	return refseq


def call_vep_and_load_info(vcfFile, keyword="none", cutoff='',
				doCallVep=True, returnVepByLocus=True,
				vep_filename=False,
				):
        """Calls VEP on a vcf file and returns a table of its output.
        vcfFile -- the full path to a vcf file
        System calls VEP, opens the output vep.tmp file. Creates
        a table vepByLocus with key=locus and value=vep line. Because
        there can be multiple variants at a locus, the lookup value of a
        locus here does not give pathogenicity information.
        Returns the table of vepByLocus information.
        """
	if(keyword == "none"):
		keyword = vcfFile
	frame = inspect.currentframe()
	print "in function to call vep with %s" % str(inspect.getargvalues(frame))
	if(doCallVep):
		splitLen = 200 # 500 lines per file
		if(not os.path.exists('vcf/')):
			os.system('mkdir vcf')
		outputBase = './vcf/splitvcf.'
		initialF = open(vcfFile, 'r')
		count = 0
		at = 0
		splitVcfF = None
		splitVcfsList = set()
		for line in initialF:
		    if count % splitLen == 0:
			if splitVcfF: splitVcfF.close()
			at += 1
			splitVcfF = open(outputBase + str(at) + '.vcf', 'w')
			splitVcfsList.add(outputBase + str(at) + '.vcf')
		    splitVcfF.write(line)
		    count += 1
		initialF.close()
		splitVcfF.close()
		for splitVcfF in splitVcfsList:
			cmdl = ["perl", path_to_vep,
				"--database", "-i", splitVcfF, "--force_overwrite", 
				"-o", "%s.vep" % splitVcfF, "--symbol",
				"canoncical", "--sift", "b"]
			print str(cmdl)
			try:
				vepSub.kill()
			except:
				pass
			vepSub = subprocess.Popen(cmdl, stdout=subprocess.PIPE)
			vepSub.wait()
		concatCmd = []
		for splitVcfF in splitVcfsList:
			concatCmd.append("%s.vep" % splitVcfF)
		#fullvep = "./hets/%s.vep" % keyword
		if(vep_filename):
			fullvepf = vep_filename
		if((keyword == 'none') or (keyword == 'all')):
			fullvepf = "hets/%s/all.%s.vep" % (cutoff, cutoff)
		concatCmd = ['cat'] + concatCmd 
		print concatCmd
		with open(fullvepf, "w") as outfile:
			catSub = subprocess.Popen(concatCmd, stdout=outfile)
			catSub.wait()
	if(returnVepByLocus):
		if(keyword == 'all'):
			vep = open("./hets/%s/all.%s.vep" % (cutoff, cutoff), 'r')
		else:
			vep = open("%s.fullvep" % keyword, 'r')
		vepByLocus = dict()
		for v in vep:
			v = v.rstrip('\n')
			if(len(v.split('\t')) > 1):
				m = re.search(r'MT:(\d+)\t(\w)\t', v)
				if(m is not None):
					vepByLocus.setdefault( int(m.groups()[0]), [v] ).append( v )
		vep.close()
		for locus in vepByLocus:
			vepByLocus[locus] = set(vepByLocus[locus])
		return vepByLocus
	else:
		if(keyword == 'all'):
			return "hets/%s/all.%s.vep" % (cutoff, cutoff)
		return "%s.fullvep" % keyword


def update_with_vep_output(fileToUpdate, vcfFilename, updatedFile, keyword='none', doCallVep=True):
        """Calls VEP and updates an input file with VEP information.
        fileToUpade -- an input file with the first line as a locus
        vcfFilename -- a vcf file with the alleles covered in the file to update
        updatedFile -- filename of the file to create with vep information added
        Initially calls callVepOnVcf() and reads in the table of vep information.
        The table of vep information is used to lookup the locus at each line
        of the fileToUpdate. VEP lines containing 'upstream'/'downstream' 
        are ignored. This function merely tags the input lines with vep lines,
        but since it does not check if the variant is the same, it does not
	check if the allele is actually pathogenic.
        """
        hetsByPerson = open(fileToUpdate, 'r')
        hetsByPersonWithVep = open(updatedFile, 'w')
	# Get a table of vep information.
        vepByLocus = call_vep_and_load_info(vcfFilename, keyword, doCallVep)
        for l in hetsByPerson:
                s = l.rstrip('\n').split('\t')
                lineOut = l.rstrip('\n')
                if(int(s[1]) in vepByLocus):
                        goodLineExists = 0
			vepLinesToOutput = set()
                        for vepLine in vepByLocus[ int(s[1]) ]:
                                if(not re.search(r'stream', vepLine)):
					vepLinesToOutput.add(vepLine)
                                        goodLineExists = 1
			if(len(vepLinesToOutput)>0):
				for vepLine in vepLinesToOutput:
                                        lineOut += "\t$" + vepLine
                lineOut += "\n"
                hetsByPersonWithVep.write(lineOut)
        hetsByPerson.close()
        hetsByPersonWithVep.close()


def identifyInterestingVepLines(vepfile, outvepfile):
	vepf = open(vepfile, 'r')
	vepfout = open(outvepfile, 'w')
	for li in vepf:
		s = li.rstrip('\n').split('\t')
		_stream = re.search(r'upstream_gene_variant|downstream_gene_variant', li)
		_tol = re.search(r'tolerated', li)
		_syn = re.search(r'synonymous', li)
		_noncoding = re.search(r'non_coding', li)
		if(_stream is None and _tol is None and 
			_syn is None and _noncoding is None):
			vepfout.write(li)
	vepf.close()
	vepfout.close()


def readInterestingVepLines(vepfilename):
	vepf = open(vepfilename, 'r')
	vep_by_variant = dict()
	for li in vepf:
		if(re.match(r'\A#', li)):
			continue
		m = re.search(r'MT:(\d+)\t(\w)\t', li)
		if(m is not None):
			locus = int(m.group(1))
			mutant_variant = m.group(2)
			vep_by_variant.setdefault((locus, mutant_variant), [li]).append(li)
	vepf.close()
	for allele in vep_by_variant:
		vep_by_variant[allele] = list(set(vep_by_variant[allele]))[0]
	return vep_by_variant


def writeVcf(hets_filename, justReturnFilename, refseq):
	"""Write .vcf from hets file.
	hets_filename -- file of heteroplasmies containing the locus 
	on column two and a bc string somewhere in the line.
	Returns the .vcf filename.
	"""
	shortPathMatch = re.match(r'\A([^/]+)$', hets_filename)
	if(shortPathMatch is not None):
		outBasename = shortPathMatch.group(1)
	else:
   		outBasename = re.match(r'.*/([^/]+)$', hets_filename).group(1)
	if(justReturnFilename):
		return "%s.vcf" % outBasename
	allVcf = open("%s.vcf" % outBasename, 'w')
        majorVcf = open("%s.hetMajor.vcf" % outBasename, 'w')
        minorVcf = open("%s.hetMinor.vcf" % outBasename, 'w')
	hetsF = open(hets_filename, 'r')

	for li in hetsF:
		s = li.rstrip('\n').split('\t')
		locus = s[1]
		bcre = re.search(r'\[\(\'(\w)\', (\d+)\), \(\'(\w)\', (\d+)\), \(\'(\w)\', (\d+)\), \(\'(\w)\', (\d+)\)', li).groups()
		bc = { bcre[0]: bcre[1], bcre[2]: bcre[3], bcre[4]: bcre[5], bcre[6]: bcre[7] }
		sorted_bc = sorted(bc.items(),key=itemgetter(1))

		if( (sorted_bc[-1][0] != refseq[int(locus)-1])
		    and (sorted_bc[-1][1] > 0)):
			vcfLine = "MT\t%i\t.\t%s\t%s\t999\t.\tDP=1000\n" % (
			int(locus), refseq[int(locus)-1] , sorted_bc[-1][0] )
			majorVcf.write(vcfLine)
			allVcf.write(vcfLine)
		if( (sorted_bc[-2][0] != refseq[int(locus)-1])
		    and (sorted_bc[-2][1] > 0)):
			vcfLine = "MT\t%i\t.\t%s\t%s\t999\t.\tDP=1000\n" % (
			int(locus), refseq[int(locus)-1] , sorted_bc[-2][0] )
			minorVcf.write(vcfLine)
			allVcf.write(vcfLine)
        allVcf.close()
        majorVcf.close()
        minorVcf.close()

	return "%s.vcf" % outBasename


def getGenbankFrequencies(gb_ref_filename='./lib/mitobank_variant_freqs.txt'):
	"""Reutrn a dict of the genbank frequency of all alleles.
	The dict returned holds key=(int(locus), mutant_variant) tuples
	"""
	src_path = os.path.dirname(os.path.realpath(__file__))
	gb_ref_filename = "%s/lib/mitobank_variant_freqs.txt" % src_path
	gb_frequencies = dict() # key=(locus, variant) tuple
	if(os.path.exists(gb_ref_filename)):
		with  open(gb_ref_filename, 'r') as gb_ref:
			next(gb_ref)
			for li in gb_ref:
				s = li.rstrip('\n').split('\t')
				allele = (int(s[0]), s[2])
				gb_frequencies[allele] = int(s[3])
		return gb_frequencies
	else:
		return False


def readMitobank(filename='lib/MutationsCodingControl.mitomap.txt'):
	"""Reads information from mitomap files and reutrns dict.
	The dict returned holds key=(int(locus), mutant_variant) tuples
	"""
	src_path = os.path.dirname(os.path.realpath(__file__))
	filename = "%s/lib/MutationsCodingControl.mitomap.txt" % src_path
        variants = dict()
        mitobankf = open(filename, 'r')
        for li in mitobankf:
                s = li.rstrip('\n').split('\t')
                try:
                        locus = s[0]
                        m = re.match('([ACGTB])(\d+)([ACGTN])$', s[3])
                        if(m is not None):
                                wt = m.group(1)
                                mut = m.group(3)
                                allele = (int(locus), mut)
                                variants[allele] = li.rstrip('\n')
                except:
                        print "error parsing line" + li
        mitobankf.close()
	filename = "%s/lib/MutationsRNA.txt" % src_path
        mitobankf = open(filename, 'r')
        for li in mitobankf:
                s = li.rstrip('\n').split('\t')
                try:
                        locus = s[0]
                        m = re.match('([ACGTB])(\d+)([ACGTN])$', s[3])
                        if(m is not None):
                                wt = m.group(1)
                                mut = m.group(3)
                                allele = (int(locus), mut)
                                variants[allele] = li.rstrip('\n')
                except:
                        print "error parsing line" + li
        mitobankf.close()
        return variants


def read_missense(filename=False):
	"""Reads information from mitomap files and reutrns dict.
	The dict returned holds key=(int(locus), mutant_variant) tuples
	"""
	if not filename:
		src_path = os.path.dirname(os.path.realpath(__file__))
		filename = src_path + '/lib/pathogenic_alleles/missense.txt'
		print "missense path = %s" % filename
        variants = dict()
        miss_f = open(filename, 'r')
	next(miss_f)
        for li in miss_f:
                s = li.rstrip('\n').rstrip('\r').split('\t')
                try:
			# s[1] = A>C format
                        locus = s[0]
                        m = re.match('([ACGT])>([ACGT])', s[1] )
                        if(m is not None):
                                wt = m.group(1)
                                mut = m.group(2)
                                allele = (int(locus), mut)
                                #variants[allele] = li.rstrip('\n')
				variants[allele] = {'mutpred_score': s[8],
						'mtDNA_selection_score': s[9]}
                except:
                        print "error parsing line" + li
	miss_f.close()
	return variants


def read_trna(filename=False):
	if not filename:
		src_path = os.path.dirname(os.path.realpath(__file__))
		filename = src_path + '/lib/pathogenic_alleles/all_trna.txt'
		print "missense path = %s" % filename
        trna_f = open(filename, 'r')
	next(trna_f)
	variants = {}
        for li in trna_f:
                s = li.rstrip('\n').split('\t')
                try:
                        locus = s[3]
			mut = s[2]
                        allele = (int(locus), mut)
                        variants[allele] = li.rstrip('\n')
                except:
                        print "error parsing line" + li
        trna_f.close()
        return variants

def load_omim(path_to_omim='lib/OMIMpathogenic'):
        # Check OMIM for information on pathogenic variants.
        omim = dict()
        if(os.path.exists(path_to_omim)):
                omimF = open(path_to_omim, 'r')
                for l in omimF:
                        s = l.strip('\n').split('\t')
                        omim[ str(s[0]) ] = { 'pathogenic': s[1], 'info': str(s[2:]) }
                omimF.close()
                return omim
        else:
                omim = False


def is_known_heteroplasmy(loci, file_id, keyword):
	# Can we find the right hets file?
	loci = int(loci)
	if(not os.path.exists("%s.fileids.hets" % keyword)):
		print "hettools.is_known_heteroplasmy:"
		print "%s.fileids.hets file is missing!" % keyword
		return False
	hetsf = open("%s.fileids.hets" % keyword, 'r')
	for li in hetsf:
		s = li.rstrip('\n').split('\t')
		if((s[0] == file_id) and (int(s[1])==loci)):
			return True
	return False
	hetsf.close()


def read_phylotree(filename=r'lib/phylotreeAlleles.txt'):
	"""
	phylo dict key=(int locus, variant) tuple format. Values are all 1.
	"""
	src_path = os.path.dirname(os.path.realpath(__file__))
	filename = "%s/lib/phylotreeAlleles.txt" % src_path
	if(not os.path.exists(filename)):
		print "Phylotree file missing from %s." % filename
		return False
	phylo = dict()
	with open(filename, 'r') as f:
		for li in f:
			m = re.match(r'\A(\d+)([a-zA-Z]+)$', li)
			if m is not None:
				phylo[(int(m.group(1)), m.group(2))] = 1
			else:
				print "Unexpected format in phylotree line: %s" % li
	return phylo


def get_all_individuals(paths_fname=r'./lib/paths2.paths', get_id_to_filename=False):
	het_p = set()
	het_m = set()
	het_f = set()
	paths = load_paths(paths_fname, 'all', getAll=True)
	idToFilename = {}
	bcFilenameToId = {}
	probands = dict()
	mothers = dict()
	fathers = dict()
	families = dict()
	for path in paths:
		(_idToFilename, _bcFilenameToId) = read_map_file(path['map_file'])
		(_probands, _mothers, _fathers, _families, _unused_ids
			) = read_in_stableid_families(path['trios_file'])
		idToFilename.update(_idToFilename)
		bcFilenameToId.update(_bcFilenameToId)
		probands.update(_probands)
		mothers.update(_mothers)
		fathers.update(_fathers)
		families.update(_families)
	probands = set(probands)
	mothers = set(mothers)
	fathers = set(fathers)
	if(get_id_to_filename):
		return [probands, mothers, fathers, bcFilenameToId, idToFilename]
	return [probands, mothers, fathers, bcFilenameToId]


def concat_hets(keyword_list=['init', 'y2', 'y3'], force=False):
	for cutoff in ['high', 'mid', 'low']:
		fname = "hets/%s/all.%s.concordance.all" % (cutoff, cutoff)
		cat_het_fname = "hets/%s/all.%s.hets" % (cutoff, cutoff)
		if((not os.path.exists(cat_het_fname))
			or force):
			cmdl = ['cat']
			for keyword in keyword_list:
				if(cutoff == 'high'):
					a_het_fname = "hets/%s/%s.hets" % (cutoff, keyword)
				else:
					a_het_fname = "hets/%s/%s.%s.hets" % (
					cutoff, keyword, cutoff)
				cmdl += [a_het_fname]
			with open(cat_het_fname, 'w') as cat_het:
				print cmdl
				subP = subprocess.Popen(cmdl, stdout=cat_het)
				subP.wait()
	return True
