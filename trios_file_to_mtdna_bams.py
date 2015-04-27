"""
Process a file of trios.
Splits by folder name, based on the proband's .bam folder location.
Extracts the mtDNA, checks for complete trios, writes split .trios, .map and .path files.
Ouputs the bam_to_hets.py command needed to process the mtNDA.

Example:
	python trios_file_to_mtdna_bams.py trios_file

Trios files are in the format:
proband_id\tmother_id\tfather_id\tproband_path\tmother_path\tfather_path\n
Use something like:
awk -F$'\t' 'BEGIN {OFS=FS} {print $1, $2, $3, $10, $11, $12}' triosfile > formatted_triosfile
to reformat trios files into the required format.

"""
import re
import sys
import os
import subprocess


class trio():
	def __init__(self, cols):
		self.cols = cols
		# 0-based col 9 is the proband's bam file path.
		self.proband_path = os.path.dirname(self.cols[3].rstrip('\r'))
		self.output = os.path.basename(os.path.dirname(self.proband_path))
		self.paths = [
			os.path.dirname(self.cols[3]),
			os.path.dirname(self.cols[4]),
			os.path.dirname(self.cols[5])]

	def print_line(self):
		return "\t".join(self.cols) + "\n"

	def write_map_lines(self):
		li = "%s\t%s\n" % (self.cols[0], self.to_bc_filename(self.cols[3]))  # Proband.
		li += "%s\t%s\n" % (self.cols[1], self.to_bc_filename(self.cols[4]))
		li += "%s\t%s\n" % (self.cols[2], self.to_bc_filename(self.cols[5]))
		return li

	def to_bc_filename(self, bam_filename):
		return os.path.basename(bam_filename).rstrip(r'.bam') + '.mt.bam.bc'

def extract_mtdna(input_bamfile, output_directory):
	if not os.path.exists(input_bamfile):
		#print "Does not exist: %s" % input_bamfile
		return False
	output_directory = output_directory.rstrip('\n') + '/'
	if not os.path.exists(output_directory):
		os.system('mkdir %s' % output_directory)
		return True
	output_name = output_directory + os.path.basename(input_bamfile).replace('\r', '')
	output_name = output_name.rstrip(r'.bam') + '.mt.bam'
	print "output filename = %s" % output_name
	cmdl = ['samtools', 'view', '-h', '-b', input_bamfile,
		'-q', '20', 'MT']
	print cmdl
	with open(output_name, 'w') as f:
		sub = subprocess.Popen(cmdl, stdout=f)
		sub.wait()
	cmdl = ['samtools', 'sort', output_name, output_name.rstrip('.bam')]
	print cmdl
	os.system(" ".join(cmdl))
	cmdl = ['samtools', 'index', output_name]
	os.system(" ".join(cmdl))

def create_map_file(trios_at_path, out_folder):
	output_map_filename = 'het_burden/lib/%s.map2' % out_folder
	with open(output_map_filename, 'w') as f:
		for _trio in trios_at_path:
			f.write(_trio.write_map_lines())

def create_paths_file(trios_at_path, out_folder):
	output_trios_filename = 'het_burden/lib/%s.paths' % out_folder
	with open(output_trios_filename, 'w') as f:
		f.write("dummy paths line")

def create_trios_file(trios_at_path, out_folder):
	output_trios_filename = 'het_burden/lib/%s.trios2' % re.sub('/', '', out_folder)
	print "output trios file: %s" % output_trios_filename
	output_trios_file = open(output_trios_filename, 'w')
	for _trio in trios_at_path:
		output_trios_file.write(_trio.print_line())
	output_trios_file.close()

def count_num_that_exist(trios_at_path, file_exists, _path):
	for _trio in trios_at_path:
		file_exists[_path][_trio] = True
		for bamfile in _trio.cols[3:]:
			if not os.path.exists(bamfile):
				file_exists[_path][_trio] = False
	without_file = len([x for x in file_exists[_path] if not file_exists[_path][x]])
	with_file = len([x for x in file_exists[_path] if file_exists[_path][x]])
	return (with_file, without_file)

def extract_mtdna_from_all_trios(trios_at_path, out_folder):
	for _trio in trios_at_path:
		for bamfile in _trio.cols[3:]:
			extract_mtdna(bamfile, out_folder)

trios = []
trios_filename = sys.argv[1]  #'./het_burden/lib/all_hgi_improved_bams_for_4344_trios.txt'
with open(trios_filename, 'r') as f:
	next(f)
	for li in f:
		s = li.rstrip('\n').rstrip('\r').split('\t')
		trios.append(trio(s))

paths = set()
paths_n = {}
for _trio in trios:
	# 9 to 11, the last three, are full paths.
	paths.add(os.path.dirname(_trio.cols[3]))
	paths.add(os.path.dirname(_trio.cols[4]))
	paths.add(os.path.dirname(_trio.cols[5]))

for _path in paths:
	paths_n[_path] = 0

for _trio in trios:
	paths_n[os.path.dirname(_trio.cols[3])] += 1
	paths_n[os.path.dirname(_trio.cols[4])] += 1
	paths_n[os.path.dirname(_trio.cols[5])] += 1

print "Trios: %i Trios * 3: %i" % (len(trios), int(3 * len(trios)))
print "Found paths: "
for _path in paths:
	print _path

print "Counts: "
for _path in paths_n:
	print "%s: %i" % (_path, paths_n[_path])
print "Sum of counts: %i" % (sum(
	[int(paths_n[_path]) for _path in paths_n.keys()]
	))
# Output folder is named for the proband folder.
file_exists = {}
do_process = set()
for _path in paths:
	file_exists[_path] = {}
	trios_at_path = [_trio for _trio in trios if _trio.proband_path == _path]
	print "%s: Trios: %i." % (_path, len(trios_at_path))
	in_folder = _path
	out_folder = trios_at_path[0].output
	(with_file, without_file) = count_num_that_exist(trios_at_path, file_exists, _path)
	print "%s: With files: %i. Without files: %i. Total %i." % (
		_path, with_file, without_file, len(file_exists[_path]))
	#extract_mtdna_from_all_trios(trios_at_path, out_folder)
	if with_file > 0: #and without_file == 0:
		do_process.add((_path, out_folder))
		create_trios_file(trios_at_path, out_folder)
		create_map_file(trios_at_path, out_folder)
		create_paths_file(trios_at_path, out_folder)
	
print "Run the following command to process the mtDNA:"
li = "python bams_to_hets.py"
for (_path, keyword) in do_process:
	li += " %s" % keyword
print li


