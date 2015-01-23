"""
Calls heteroplasmies on a list of bam folders.

Example: python bams_to_hets.py bamfolder1 bamfolder2

Expects:
A het_burden/lib/bamfoldername.path file for each bam folder.
A het_burden/lib/bamfoldername.trios2 file for each bam folder.
A het_burden/lib/bamfoldername.map2 file for each bam folder.

"""
import os
import re
import sys
import argparse

if __name__ == "__main__":
	bam_folders = sys.argv[1:len(sys.argv)]
	bc_folders = set()
	for bam_folder in bam_folders:
		bc_folder = re.sub("bams?", "", bam_folder)
		bc_folder += "bc" 
		cmd = "python findHets/write_basecall_files.py "
		cmd += "-i %s -o %s" % (bam_folder, bc_folder)
		print cmd
		os.system(cmd)
		bc_folders.add(bc_folder)
	for bc_folder in bc_folders:
		cmd = "python findHets/findHets.py "
		cmd += "-b %s -o %s" % (bc_folder, bam_folder)
		print cmd
		os.system(cmd)
	for bamfolder in bam_folders:
		bc_folder = re.sub("bams?", "", bam_folder)
		bc_folder += "bc" 
		cmd = "python het_burden/het_burden.py "
		cmd += "-k %s -p het_burden/lib/%s.path" % (
			bc_folder, bam_folder)
		print cmd
		os.system(cmd)
		cmd = "python het_burden/het_concordance_tools.py "
		cmd += "-p het_burden/lib/%s.path -k %s"
		print cmd
		os.system(cmd)
		cmd = "python het_burden/categorize_alleles.py "
		cmd += "-d ./hets"
		print cmd
		os.system(cmd)
