Heteroplasmy Inheritance
==================

Identifies heteroplasmy in family trios - a mother, a father, and a single offspring.

extract_unique_mt.sh
---------
extract_unique_mt.sh calls samtools to extract MT reads.

A -q 20 MAPQ cutoff is applied, which serves to select for unique sequences.

Example: ::
	$ bash extract_unique_mt.sh infolder outfolder


write_basecall_files.py
---------
write_basecall_files converts a folder of .bam files into a folder of .bc files.

.bc files represent reads at each locus of each base. This code is intended to replicate Mitoseek in Python.

Example: ::

	$ python write_basecall_files -i bamfolder -o bcfolder -f lib/rcrs.fasta

findHets.py
-----------
findHets.py finds heteroplasmic individuals in a folder of .bc files.

Example: ::

	$ python findHets.py -b bcfolder -o prefix_for_output_files.

This program will put results in a ./hets folder, separated by cutoff.
