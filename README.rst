Heteroplasmy Inheritance
==================

Identifies and analyzes heteroplasmy in family trios - a mother, a father, and a single offspring.


extract_unique_mt.sh
---------
extract_unique_mt.sh calls samtools to extract MT reads.

A -q 20 MAPQ cutoff is applied, which serves to select for unique sequences.

Example: ::

	$ bash extract_unique_mt.sh example_1000_genomes_bams/ outfolder

bams_to_hets.py
---------
Identifies heteroplasmy and analyzes alleles for a folder of bam files.

All scripts listed below this script are included within it.

This script is a complete pathway, which expects certain library files to exist.

Other scripts in this package do not require such library files and can be run on their own.

Example: ::

	$ python bams_to_hets.py bamfolder1 bamfolder2 bamfolder3

Expects the following files to exist for each bam folder:

het_burden/lib/bamfolder.path
het_burden/lib/bamfolder.map2
het_burden/lib/bamfolder.trios2

The format of these folders is given by example in het_burden/lib/1kgenomes.*

Briefly, .path links files and folders, .map2 relates individuals to files, and .trios2 describes family relations.

See formats.txt for specifications.

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


het_burden
---------
Scripts in het_burden produce analytical information for ./hets folders.

These scripts require lib/.path, lib/.map2 and lib/.trios2, specifying relations between individuals.

These scripts also require allele annotation files (obtained from Mitobank) included in lib/. tRNA mutation information is taken from "Prediction of pathogenic mutations in mitochondrially encoded human tRNAs." (F. Kondrashov, 2005).
 
Example: ::

	$ python het_burden/het_burden.py -k 1kgenomes -p het_burden/lib/1kgenomes.path

...will output various metrics to ./hets/burden/cutoff.

Example: ::

	$ python het_burden/het_concordance_tools.py -p het_burden/lib/1kgenomes.path -k 1kgenomes

...will output ./hets/burden/cutoff/keyword.concordance files containing information on the concordance of alleles between mother and offspring.

het_burden scripts are a mess.
