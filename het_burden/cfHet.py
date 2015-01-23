import re
import os
from operator import itemgetter
import sys


class cfHet:
	"""Holds information on heteroplasmies for a given mother-proband pair.
	self.name = (proband, mother) tuple that identifies the instance
	self.loci = dict of key=loci; value=dict of mother/proband;
	loci[locus][p/m] = tuple of (base observations, minor allele frequency)
	"""
	def __init__(self, proband="", mother="", father=""):
		self.proband = proband
		self.mother = mother
		self.father = father
		self.name = (proband, mother)
		# information on heterozygous loci. All dict objects are key=loci.
		self.loci = dict()
		self.concord = dict()
		self.alleles = dict() #not used
		# information from vep. key=loci.
		self.siftType = dict()
		self.siftLine = dict()
		# information for outputing for statistics
		self.frequencyTuples = dict()
		self.frequency_tuples_alleles = dict() 
		# ^ holds what freq is what allele to compare vs the sift allele

	def fill_in_family(self, stable_id, lookup_family, families):
		if(self.proband == ""):
			self.proband = families[lookup_family[stable_id]][0]
		if(self.mother == ""):
			self.mother = families[lookup_family[stable_id]][1]
		if(self.father == ""):
			self.father = families[lookup_family[stable_id]][2]

	def parse_sift(self, locus, base_observations, vep):
		vep = re.sub('\n', '', vep)
		# vep is an empty string if there is no vep information in the file
		if(vep):
			vepSplit = vep.rstrip('\n').split('\t')
			allele = (locus, vepSplit[2])
			sift_regex = re.search(r'SIFT=(\w+)\(([0-9\.]+)\)', vep)
			if(sift_regex is not None):
				self.siftType[allele] = sift_regex.groups()[0]
				self.siftLine[allele] = vep.rstrip('\n')

	def add_depths(self, bcdir, idToFilename):
		self.depth = dict()
		if(self.proband in idToFilename):
			bcfile = bcdir.rstrip('/') + '/' + idToFilename[self.proband]
			self.depth['proband'] = self.add_depth(
				self.proband, bcfile)
		if(self.mother in idToFilename):
			bcfile = bcdir.rstrip('/') + '/' + idToFilename[self.mother]
			self.depth['mother'] = self.add_depth(
				self.mother, bcfile)
		if(self.father in idToFilename):
			bcfile = bcdir.rstrip('/') + '/' + idToFilename[self.father]
			self.depth['father'] = self.add_depth(
				self.father, bcfile)

	def add_depth(self, stable_id, bcfile):
		totaldepth = 0
		if(os.path.exists(bcfile)):
			with open(bcfile, 'r') as bcfileh:
				for li in bcfileh:
					# chr locus ref_allele basecounts
					q = li.rstrip('\n').split('\t')
					for i in q[3:11]:
						totaldepth += int(i)
		else:
			print "no such file %s" % bcfile
		totaldepth = float(totaldepth)/float(16569)
		return totaldepth

	def set_locus_from_hets_file(self, locus, member, maf, base_observations, vep):
		"""Takes the arguments:
		locus, stable_id of memeber, minor allele frequency, dict of observations, vep
		Checks if the stable_id is the mother or proband of this object instance.
		Sets the self.loci dict at the given locus, under self.loci[locus][p or m]
		to a tuple of maf and the dict of base_observations.
		Sets the self.concord[locus][p/m maf/evidence] objects.
		The self.loci and self.concord objects are used to establish concordance.
		""" 
		locus = int(locus)
		maf = float(maf)
		self.parse_sift(locus, base_observations, vep)
		for k in base_observations:
			base_observations[k] = int(base_observations[k])
		if(member == self.proband):
			if(locus in self.loci):
				self.loci[locus]['proband'] = (maf, base_observations) 
			else:
				self.loci[locus] = {'proband': (maf, base_observations)}
			self.concord.setdefault(locus,
				 {'proband_maf': float(maf)})['proband_maf'] = float(maf)
			self.concord.setdefault(locus, 
				{'proband_evidence': 'strong'})['proband_evidence'] = 'strong'
		if(member == self.mother):
			if(locus in self.loci):
				self.loci[locus]['mother'] = (maf, base_observations) 
			else:
				self.loci[locus] = {'mother': (maf, base_observations)}
			self.concord.setdefault(locus, 
				{'mother_maf': float(maf)})['mother_maf'] = float(maf)
			self.concord.setdefault(locus, 
				{'mother_evidence': 'strong'})['mother_evidence'] = 'strong'
		if(member == self.father):
			if(locus in self.loci):
				self.loci[locus]['father'] = (maf, base_observations) 
			else:
				self.loci[locus] = {'father': (maf, base_observations)}

	def set_locus_from_bc_file(self, locus, member, maf, evidence, base_observations):
		"""Takes the arguments:
		locus, stable_id of memeber, minor allele frequency, evidence,  dict of observations
		Checks if the stable_id is the mother or proband of this object instance.
		Sets the self.loci dict at the given locus, under self.loci[locus][p or m]
		to a tuple of maf and the dict of base_observations.
		Sets the self.concord[locus][p/m maf/evidence] objects.
		The self.loci and self.concord objects are used to establish concordance.
		""" 
		locus = int(locus)
		maf = float(maf)
		for k in base_observations:
			base_observations[k] = int(base_observations[k])
		if(locus not in self.loci):
			self.loci[locus] = dict()
		if(locus not in self.concord):
			self.concord[locus] = dict()
		if(member == self.proband):
			self.loci[locus]['proband'] = (maf, base_observations)
			self.concord[locus]['proband_maf'] = maf
			self.concord[locus]['proband_evidence'] = evidence
		if(member == self.mother):
			self.loci[locus]['mother'] = (maf, base_observations)
			self.concord[locus]['mother_maf'] = maf
			self.concord[locus]['mother_evidence'] = evidence
		if(member == self.father):
			self.loci[locus]['father'] = (maf, base_observations)

	def is_above_depth(self,locus):
		locus = int(locus)
		mother_dictOfBaseObservations = self.loci[locus]['mother'][1]
		proband_dictOfBaseObservations = self.loci[locus]['proband'][1]
		mbc = sorted(mother_dictOfBaseObservations.items(), key=itemgetter(1))
		pbc = sorted(proband_dictOfBaseObservations.items(), key=itemgetter(1))
		mdepth = float(0)
		for i in mbc:
			mdepth += float(i[1])
		pdepth = float(0)
		for i in pbc:
			pdepth += float(i[1])
		cutoff = 10
		if((mdepth < cutoff) and (pdepth >= cutoff)): #if mother not above cutoff
			self.frequencyTuples[locus] = [('', ''), (float(pbc[-1][1])/pdepth, 
								float(pbc[-2][1])/pdepth)]
			return False
		if((pdepth < cutoff ) and (mdepth >= cutoff)): #if mother not above cutoff
			self.frequencyTuples[locus] = [(float(mbc[-1][1])/mdepth, 
							float(mbc[-2][1])/mdepth), ('', '')]
			return False
		if((mdepth < cutoff) and (pdepth < cutoff )):
			self.frequencyTuples[locus] = [('',''), ('', '')]
			return False
		return True

	def set_frequencies(self, locus):
		"""Sets the self.frequencyTuples at a given locus so that
		.concordance files can be written.
		Takes a locus as an argument. Retrieves the base observations
		from self.loci and sets self.frequencyTuples based on that 
		info.
		"""
		locus = int(locus)
		if(locus not in self.loci):
			return False
		if('mother' not in self.loci[locus]):
			return False
		if('proband' not in self.loci[locus]):
			return False
		mother_dictOfBaseObservations = self.loci[locus]['mother'][1]
		proband_dictOfBaseObservations = self.loci[locus]['proband'][1]
		mbc = sorted(mother_dictOfBaseObservations.items(), key=itemgetter(1))
		pbc = sorted(proband_dictOfBaseObservations.items(), key=itemgetter(1))
		mdepth = float(0)
		for i in mbc:
			mdepth += float(i[1])
		pdepth = float(0)
		for i in pbc:
			pdepth += float(i[1])
		cutoff = 10
		if((mdepth < cutoff) and (pdepth >= cutoff)):  # if mother below cutoff
			self.frequencyTuples[locus] = [('', ''), (float(pbc[-1][1])/pdepth, 
								float(pbc[-2][1])/pdepth)]
		if((pdepth < cutoff ) and (mdepth >= cutoff)):  # if proband below cutoff
			self.frequencyTuples[locus] = [(float(mbc[-1][1])/mdepth, 
							float(mbc[-2][1])/mdepth), ('', '')]
		if((mdepth < cutoff) and (pdepth < cutoff)):  # both below cutoff
			self.frequencyTuples[locus] = [('',''), ('', '')]
		# is this a simple case of there being only two alleles in play?
		allelesInPlay = set()
		for j in [mbc[-1], mbc[-2], pbc[-1], pbc[-2]]:
			if(int(j[1]) > 1):
				allelesInPlay.add(j[0])
		# avoid division by zero errors in the next section
		# TO DO: a better solution
		if(mdepth==0):
			mdepth=1
		if(pdepth==0):
			pdepth=1
		if(len(allelesInPlay) < 3):
			# Only two alleles at issue. 
			# Who has the maximum frequency of a minor allele?
			if(mbc[-2][1] >= pbc[-2][1]): 
				# mbc[-2][0] is the minor allele of interest (holds nucleotide)
				m_majorFreq = float(mbc[-1][1])/mdepth
				m_minorFreq = float(mbc[-2][1])/mdepth
				p_majorFreq = float(
					proband_dictOfBaseObservations[mbc[-1][0]])/pdepth
				p_minorFreq = float(
					proband_dictOfBaseObservations[mbc[-2][0]])/pdepth
				# The mother has the greater minor allele frequency.
				# We set the alleles tuple to the (major, minor) allele
				# of the mother in this case. The mothers major/minor alleles
				# remain as-is in the frequencyTuple object. The proband's 
				# frequencies are given in the same order, even if the "minor"
				# proband allele is now the more abundant allele.
				self.alleles[locus] =  (mbc[-1][0], mbc[-2][0])
				self.frequency_tuples_alleles[locus] = (mbc[-1][0], mbc[-2][0]) 
			if(mbc[-2][1] < pbc[-2][1]):
				m_majorFreq = float(
					mother_dictOfBaseObservations[pbc[-1][0]])/mdepth
				m_minorFreq = float(
					mother_dictOfBaseObservations[pbc[-2][0]])/mdepth
				p_majorFreq = float(pbc[-1][1])/pdepth
				p_minorFreq = float(pbc[-2][1])/pdepth
				self.alleles[locus] =  (pbc[-1][0], pbc[-2][0])
				self.frequency_tuples_alleles[locus] = (pbc[-1][0], pbc[-2][0])
			# Set based on who has the greater minor allele frequency
			self.frequencyTuples[locus] = [
				(m_majorFreq, m_minorFreq),
				(p_majorFreq, p_minorFreq)]
		else: 
		# Major allele is set to be the mother's by default. 
		# Minor allele is set to be the highest frequency minor allele of 
		# either mother or proband.
			print "tri-allelic situation: "
			print "%s. mbc: %s pbc: %s" % (str(allelesInPlay), str(mbc), str(pbc))
			if( mbc[-2][1] >= pbc[-2][1] ): # mbc[-2][0] is the minor allele of interest
				self.alleles[locus] =  (mbc[-1][0], mbc[-2][0])
				m_majorFreq = float(mbc[-1][1])/mdepth
				m_minorFreq = float(mbc[-2][1])/mdepth
				p_majorFreq = float(proband_dictOfBaseObservations[mbc[-1][0]])/pdepth
				p_minorFreq = float(proband_dictOfBaseObservations[mbc[-2][0]])/pdepth
				self.frequency_tuples_alleles[locus] = (mbc[-1][0], mbc[-2][0]) 
			if( mbc[-2][1] < pbc[-2][1] ):
				self.alleles[locus] =  (pbc[-1][0], pbc[-2][0])
				m_majorFreq = float(mbc[-1][1])/mdepth
				m_minorFreq = float(mother_dictOfBaseObservations[pbc[-2][0]])/mdepth
				p_majorFreq = float(proband_dictOfBaseObservations[mbc[-1][0]])/pdepth
				p_minorFreq = float(pbc[-2][1])/pdepth
				self.frequency_tuples_alleles[locus] = (mbc[-1][0], pbc[-2][0]) 
			self.frequencyTuples[locus] = [
				(m_majorFreq, m_minorFreq),
				(p_majorFreq, p_minorFreq)]
		return True

	def return_alleles_in_order(self, locus):
		locus = int(locus)
		print "setting freq value for the locus -%s-" % str(locus)
		if( (not hasattr(self, 'frequency_tuples_alleles')) or 
			(locus not in self.frequency_tuples_alleles)):
			self.set_frequencies(locus)
		print "returning the value -%s-" % str(self.frequency_tuples_alleles[locus])
		return self.frequency_tuples_alleles[locus]

	def write_frequencies(self, locus, pathogenicity, mode='mitobank', mitobank={}):
		locus = int(locus)
		if((re.search('Error', self.concord[locus]['proband_evidence']) is not None) 
		 or (re.search('Error', self.concord[locus]['mother_evidence']) is not None) ):
			return False
		if(locus not in self.frequencyTuples):
			return False
		if(pathogenicity == 'all'):
			return "%f\t%f\t%f\t%f" % (self.frequencyTuples[locus][0][0], 
			self.frequencyTuples[locus][0][1], 
			self.frequencyTuples[locus][1][0], 
			self.frequencyTuples[locus][1][1])
		if(pathogenicity != 'all'):
			if(mode == 'mitobank'):
				return self.write_frequencies_mitobank(locus, mitobank)
			# if this is a valid sift type we expect a score and allele set
			return self.write_frequencies_specific_pathogenicity(
					locus, pathogenicity)
		return False

	def write_frequencies_mitobank(self, locus, mitobank):
		allele_mother_major = (int(locus), self.frequency_tuples_alleles[locus][0])
		allele_mother_minor = (int(locus), self.frequency_tuples_alleles[locus][1])
		if(allele_mother_major not in mitobank
			and allele_mother_minor not in mitobank):
			return False
		if(allele_mother_major in mitobank
			and allele_mother_minor not in mitobank):
			return "%f\t%f\t%f\t%f\t%s" % (self.frequencyTuples[locus][0][1], 
				self.frequencyTuples[locus][0][0], 
				self.frequencyTuples[locus][1][1], 
				self.frequencyTuples[locus][1][0],
				mitobank[allele_mother_major])
		if(allele_mother_major not in mitobank
			and allele_mother_minor in mitobank):
			return "%f\t%f\t%f\t%f\t%s" % (self.frequencyTuples[locus][0][0], 
			self.frequencyTuples[locus][0][1], 
			self.frequencyTuples[locus][1][0], 
			self.frequencyTuples[locus][1][1],
			self.siftLine[allele_mother_minor])
		if(allele_mother_major in mitobank
			and allele_mother_minor in mitobank):
			# We'll take the minor allele as the pathogenic one.
			# Presumably this case is rare.
			return "%f\t%f\t%f\t%f\t%s" % (self.frequencyTuples[locus][0][0], 
			self.frequencyTuples[locus][0][1], 
			self.frequencyTuples[locus][1][0], 
			self.frequencyTuples[locus][1][1],
			self.siftLine[allele_mother_minor])

	def write_frequencies_specific_pathogenicity(self, locus, pathogenicity):
		# we will always output the vep locus second
		allele_mother_major = (int(locus), self.frequency_tuples_alleles[locus][0])
		allele_mother_minor = (int(locus), self.frequency_tuples_alleles[locus][1])
		if((allele_mother_major not in self.siftType) 
		and (allele_mother_minor not in self.siftType)):
		# this heteroplasmy is not of any selected type
			return False
		if((allele_mother_major in self.siftType)
		and (allele_mother_minor not in self.siftType)):
			if(self.siftType[allele_mother_major] == pathogenicity):
				# only the mother's major allele is of the selected type
				# reverse the output order so the pathogenic allele is second for mother and proband
				return "%f\t%f\t%f\t%f\t%s" % (self.frequencyTuples[locus][0][1], 
				self.frequencyTuples[locus][0][0], 
				self.frequencyTuples[locus][1][1], 
				self.frequencyTuples[locus][1][0],
				self.siftLine[allele_mother_major])
		if((allele_mother_major not in self.siftType) 
		and (allele_mother_minor in self.siftType)):
			if(self.siftType[allele_mother_minor] == pathogenicity):
				# vep allele is the mother's minor allele. Output the minor allele second like normal
				return "%f\t%f\t%f\t%f\t%s" % (self.frequencyTuples[locus][0][0], 
				self.frequencyTuples[locus][0][1], 
				self.frequencyTuples[locus][1][0], 
				self.frequencyTuples[locus][1][1],
				self.siftLine[allele_mother_minor])
		if((allele_mother_major in self.siftType) 
		and (allele_mother_minor in self.siftType)):
			if( (self.siftType[allele_mother_minor] == pathogenicity)
			and (self.siftType[allele_mother_major] == pathogenicity)):
				# both major and minor alleles are of the selected type.
				# how do we handle this case? We will output the mother's minor allele second,
				# as though only the minor allele were pathogenic
				print "Both major and minor alleles are of the selected sift type"
				return "%f\t%f\t%f\t%f\t%s***%s" % (self.frequencyTuples[locus][0][0], 
				self.frequencyTuples[locus][0][1], 
				self.frequencyTuples[locus][1][0], 
				self.frequencyTuples[locus][1][1],
				self.siftLine[allele_mother_major],
				self.siftLine[allele_mother_minor])
			if(self.siftType[allele_mother_minor] == pathogenicity):
				# vep allele is the mother's minor allele. Output the minor allele second like normal
				return "%f\t%f\t%f\t%f\t%s" % (self.frequencyTuples[locus][0][0], 
				self.frequencyTuples[locus][0][1], 
				self.frequencyTuples[locus][1][0], 
				self.frequencyTuples[locus][1][1],
				self.siftLine[allele_mother_minor])
			if(self.siftType[allele_mother_major] == pathogenicity):
				# only the mother's major allele is of the selected type
				# reverse the output order so the pathogenic allele is second for mother and proband
				return "%f\t%f\t%f\t%f\t%s" % (self.frequencyTuples[locus][0][1], 
				self.frequencyTuples[locus][0][0], 
				self.frequencyTuples[locus][1][1], 
				self.frequencyTuples[locus][1][0],
				self.siftLine[allele_mother_major])
		return False
			
	def write_concordance(self, locus):
		if(locus not in self.concord):
			return "%i\tna\tna\tna\tna" % (locus)
		if(('proband_evidence' not in self.concord[locus])
			or ('mother_evidence' not in self.concord[locus])):
			return "%i\tna\tna\tna\tna" % (locus)
		if(('proband_maf' not in self.concord[locus])
			or ('mother_maf' not in self.concord[locus])):
			return "%i\tna\tna\tna\tna" % (locus)
		return "%i\t%f\t%s\t%f\t%s " % (
			locus, self.concord[locus]['mother_maf'],
			self.concord[locus]['mother_evidence'],
			self.concord[locus]['proband_maf'],
			self.concord[locus]['proband_evidence'])

	def write_hets_format(self):
		li = ""
		for locus in self.loci:
			locus = int(locus)
			if('mother' not in self.loci[locus]):
				continue
			if('proband' not in self.loci[locus]):
				continue
			mother_dictOfBaseObservations = self.loci[locus]['mother'][1]
			proband_dictOfBaseObservations = self.loci[locus]['proband'][1]
			mbc = sorted(mother_dictOfBaseObservations.items(), key=itemgetter(1))
			pbc = sorted(proband_dictOfBaseObservations.items(), key=itemgetter(1))
			mdepth = float(0)
			for i in mbc:
				mdepth += float(i[1])
			pdepth = float(0)
			for i in pbc:
				pdepth += float(i[1])
			li += "%s" % self.proband
			li += "\t%i" % locus
			if(pdepth == 0):
				li += "\t0"
			else:
				li += "\t%f" % (float(pbc[-2][1])/pdepth)
			li += "\t%i" % int(pdepth)
			li += "\t%s\t%i" % (pbc[-1][0], pbc[-1][1])
			li += "\t%s\t%i" % (pbc[-2][0], pbc[-2][1])
			li += "\t%s" % str(pbc)
			li += "\n%s" % self.mother
			li += "\t%i" % locus
			if(mdepth == 0):
				li += "\t0"
			else:
				li += "\t%f" % (float(mbc[-2][1])/float(mdepth))
			li += "\t%i" % int(mdepth)
			li += "\t%s\t%i" % (mbc[-1][0], mbc[-1][1])
			li += "\t%s\t%i" % (mbc[-2][0], mbc[-2][1])
			li += "\t%s" % str(mbc)
			li += "\t%s" % str(mbc)
			if('father' not in self.loci[locus]):
				continue
			father_dictOfBaseObservations = self.loci[locus]['father'][1]
			fbc = sorted(father_dictOfBaseObservations.items(), key=itemgetter(1))
			fdepth = float(0)
			for i in fbc:
				fdepth += float(i[1])
			li += "\n%s" % self.father
			li += "\t%i" % locus
			if(fdepth == 0):
				li += "\t0"
			else:
				li += "\t%f" % (float(fbc[-2][1])/fdepth)
			li += "\t%i" % int(fdepth)
			li += "\t%s\t%i" % (fbc[-1][0], fbc[-1][1])
			li += "\t%s\t%i" % (fbc[-2][0], fbc[-2][1])
			li += "\t%s" % str(fbc)
			li += "\n"
		return li

	def count_pairs(self):
		het_loci = set()
		pm_pair = 0
		not_pm_pair = 0
		pm_above_cutoff = 0
		for locus in self.loci:
			locus = int(locus)
			het_loci.add(locus)
			if(('mother' not in self.loci[locus])
			or ('proband' not in self.loci[locus])):
				not_pm_pair += 1
				continue
			pm_pair += 1
			if('father' not in self.loci[locus]):
				continue
			mother_dictOfBaseObservations = self.loci[locus]['mother'][1]
			proband_dictOfBaseObservations = self.loci[locus]['proband'][1]
			mbc = sorted(mother_dictOfBaseObservations.items(), key=itemgetter(1))
			pbc = sorted(proband_dictOfBaseObservations.items(), key=itemgetter(1))
			mdepth = float(0)
			for i in mbc:
				mdepth += float(i[1])
			pdepth = float(0)
			for i in pbc:
				pdepth += float(i[1])
			if((mdepth > 10) and (pdepth > 10)):
				pm_above_cutoff += 1
		return [het_loci, pm_pair, not_pm_pair, pm_above_cutoff]
