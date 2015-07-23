# benchmarks

Various scripts for analyzing SAM and Blast alignments. Specifically,

+ compare_two_sam_files_AS.py:
	Compare the single best alignments in SAM format for two input tools

+ compute_stats_unique_alignments.py:
	Compute alignment statistics for a SAM alignment input file

+ filter_better_hits.py:
	Filter reads for which SSEARCH (or another tool) found a better alignment
	than the simulated alignment

+ graph_accuracy.py:
	Output a 2D plot for illustrating the sensitivity/selectivity of tools

+ graph_accuracy_3d.py:
	Output a 3D plot for illustrating the sensitivity/selectivity and run time
	of tools