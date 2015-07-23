#!/opt/python-2.7.3/bin/python

"""
usage  : If filtering,
          (1) Illumina reads (simulated using Mason, read origin in label)
              filter_better_hits.py illumina /path/to/simulated/illumina/reads.fasta ssearch_blast_tabular_alignments.txt output_file_name.txt
          (2) 454 reads (simulated using Mason, read origin in SAM file)
              filter_better_hits.py 454 /path/to/simulated/454/reads.sam ssearch_blast_tabular_alignments.txt output_file_name.txt
          (3) Ion Torrent reads (simulated using CureSim1.2, read origin in label)
              filter_better_hits.py ion /path/to/simulated/ion/reads.fasta ssearch_blast_tabular_alignments.txt output_file_name.txt

purpose: output read labels for which SSEARCH found better alignments on the genome than the ground-truth
author : Evguenia Kopylova (jenya.kopylov@gmail.com)
date   : Feb 09, 2015
usage  : python filter_better_hits.py ['illumina', '454' or 'ion'] [ground_truth_fp] [ssearch_alignments.m8] [output_fp]
         The 'ground_truth_fp' is a file containing the ground truth alignments. This is limited to SAM alignments
         (output by the Mason simulator) for Illumina reads, the simulated FASTA file for 454 reads (output by the Mason
          simulator) and the simulated FASTA file for Ion Torrent reads (output by CureSim). This script can be updated
         to support other ground-truth input files.
"""

import sys
from skbio.parse.sequences import parse_fasta


if __name__ == '__main__':

    technology = sys.argv[1]
    ground_truth_fp = sys.argv[2]
    blast_alignments_fp = sys.argv[3]
    output_fp = sys.argv[4]
    expected_mapping = {}
    written = False

    # parse Illumina ground-truth alignments
    if technology == "illumina":
        with open(ground_truth_fp, 'U') as ground_truth:
            for label, seq in parse_fasta(ground_truth):
                name = label.split()[0]
                contig = label.split()[1].split('=')[1]
                orig_begin = int(label.split()[4].split('=')[1])
                if name not in expected_mapping:
                    expected_mapping[name] = [contig, orig_begin]
                else:
                    raise ValueError("%s seen twice" % name)

    elif technology == "454":
        with open(ground_truth_fp, 'U') as ground_truth:
            for line in ground_truth:
                if line.startswith('@'):
                    continue
                name = line.split()[0]
                contig = line.split()[2]
                orig_begin = line.split()[3]
                if name not in expected_mapping:
                    expected_mapping[name] = [contig, orig_begin]
                else:
                    raise ValueError("%s seen twice" % name)

    elif technology == "ion":
        with open(ground_truth_fp, 'U') as ground_truth:
            for label, seq in parse_fasta(ground_truth):
                name = label
                contig = label.split('_')[0]
                orig_begin = label.split('_')[3]
                if name not in expected_mapping:
                    expected_mapping[name] = [contig, orig_begin]
                else:
                    raise ValueError("%s seen twice" % name)

    else:
        raise ValueError("unrecognized technology %s" % technology)
                

    # filter sequences
    with open(blast_alignments_fp, 'U') as blast_alignments:
        with open(output_fp, 'w') as filter_out:
            blast_hits = {}
            for line in blast_alignments:
                alignment = line.strip().split('\t')
                # check quality of new alignment with previous one
                if alignment[0] in blast_hits:
                    alignment_best = blast_hits[alignment[0]]
                    # if equal bitscore
                    if (alignment_best[10] == alignment[11] and not written):
                        filter_out.write("%s\n" % alignment[0])
                        written = True
                # new alignment
                else:
                    # first (best) alignment for read
                    if len(blast_hits) == 0:
                        blast_hits[alignment[0]] = alignment[1:]
                    else:
                        # new read, clear blast_hits dict from
                        # previous read
                        blast_hits.clear()
                        blast_hits[alignment[0]] = alignment[1:]
                        written = False
                    
                    # the best hit by SSEARCH does not equal to the ground-truth
                    if ((alignment[1] != expected_mapping[alignment[0]][0]) and (alignment[8] != expected_mapping[alignment[0]][1])):
                        filter_out.write("%s\n" % alignment[0])
                        written = True
