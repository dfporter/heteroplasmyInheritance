Heteroplasmy Inheritance
==================

Identifies heteroplasmy in family trios - a mother, a father, and a single offspring.

extract_unique_mt.sh
---------
extract_unique_mt.sh calls samtools to extract MT reads.

A -q 20 MAPQ cutoff is applied, which serves to select for unique sequences.

Example: ::

	$ bash extract_unique_mt.sh example_1000_genomes_bams/ outfolder

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
