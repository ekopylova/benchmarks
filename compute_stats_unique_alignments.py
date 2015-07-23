""" Compute alignment statistics such as:
        - average % query coverage (from CIGAR string)
        - average read length
        - average SW score
        - average number of indels per alignment
        - average number of N's in the reads
    for all alignments in SAM file.

    usage: python compute_stats_unique_alignments.py [tool_name] [alignments.sam] [reads.fasta]
"""

from skbio.parse.sequences import parse_fasta
import sys
import re

if __name__ == '__main__':
    software = sys.argv[1]
    s1_mapping_f = sys.argv[2] # alignments for all study
    reads = sys.argv[3] # reads aligned only by software
    min_score = sys.argv[4] # minimum SW score to pass E-value

    expected_mapping = {}
    expected_reads = 0

    # compute read length and number of N's in a read
    with open(reads, 'U') as reads_fp:
        for label, seq in parse_fasta(reads_fp):
            expected_reads += 1
            name = label.split()[0]
            length = len(seq)
            num_amb_nt = int(seq.count('N'))
            if name not in expected_mapping:
                expected_mapping[name] = [length, num_amb_nt]
            else:
                raise ValueError("reads: %s seen twice" % name)

    print "Done reading reads file."

    s1_mapping = {}

    # software 1
    with open(s1_mapping_f, 'U') as sam_fp:
        for line in sam_fp:
            if line.startswith('@'):
                continue
            name = line.split()[0]
            if name not in expected_mapping:
                continue
            #sw_score = line.split()[11].split(':')[2]
            # segemehl splits by tabs (label can have spaces)
            contig = line.split('\t')[2]
            if contig == "*":
                continue
            #orig_begin = int(line.split('\t')[3])
            sw_score = int(line.split('\t')[11].split(':')[2])
            cigar = line.split('\t')[5]
            pattern = re.compile('([MIDNSHPX=])')
            values = pattern.split(cigar)[:-1] # turn cigar into tuple of values
            paired = (values[n:n+2] for n in xrange(0, len(values), 2)) # pair values by twos
            align_len = 0
            num_indel = 0
            for pair in paired:
                l = int(pair[0]) # length of CIGAR event
                t = pair[1] # type of CIGAR event
                if (t == 'M'):
                    align_len += l
                elif (t == 'I'):
                    align_len += l
                    num_indel += 1
                elif (t == 'D'):
                    num_indel += 1
            if name not in s1_mapping:
                s1_mapping[name] = [contig, cigar, align_len, num_indel, sw_score]
            else:
                raise ValueError("%s: %s seen twice" % (software_1, name))
    s2_mapping = {}

    print "Done reading s1 file."

    # compute stats
    total_reads = len(s1_mapping)
    total_indels = 0
    total_coverage = 0
    total_amb_nt = 0
    total_sw_score = 0
    total_read_length = 0
    low_coverage = 1.0
    high_coverage = 0.0
    num_alignments_with_min_score = 0
    for read in s1_mapping:
        total_indels += s1_mapping[read][3]
        coverage = float(float(s1_mapping[read][2])/float(expected_mapping[read][0]))
        total_coverage += coverage
        if coverage < low_coverage:
            low_coverage = coverage
        elif coverage > high_coverage:
            high_coverage = coverage
        total_amb_nt += expected_mapping[read][1]
        total_sw_score += s1_mapping[read][4]
        total_read_length += expected_mapping[read][0]
        if int(s1_mapping[read][4]) >= int(min_score):
            num_alignments_with_min_score +=1

    print "Total reads: %s" % total_reads
    print "Average number of indels per alignment: %s" % float(float(total_indels)/float(total_reads))
    print "Average query coverage: %s" % float(total_coverage/float(total_reads))
    print "Low coverage = %s\tHigh coverage = %s" % (low_coverage, high_coverage)
    print "Average number of N's: %s" % float(float(total_amb_nt)/float(total_reads))
    print "Average SW score: %s" % float(float(total_sw_score)/float(total_reads))
    print "Average read length: %s" % float(float(total_read_length)/float(total_reads))
    print "Number of alignments with >= min_score (%s) = %s" % (min_score, num_alignments_with_min_score)
