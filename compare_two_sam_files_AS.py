"""Compare alignments between two SAM files (single best alignment per read required)
   usage: python compare_two_sam_files_AS.py [tool_1_name] [tool_2_name] [reads FASTA file] [tool_1_sam] [tool_2_sam] [offset]
   The final 'offset' parameter is the +-N nucleotides range within which the starting
   alignment positions must be for two tools in order to consider the alignments equal
"""

from skbio.parse.sequences import parse_fasta
import sys
import re

if __name__ == '__main__':
    software_1 = sys.argv[1]
    software_2 = sys.argv[2]
    reads = sys.argv[3]
    s1_mapping_f = sys.argv[4]
    s2_mapping_f = sys.argv[5]
    offset = int(sys.argv[6])

    expected_mapping = {}
    expected_reads = 0

    # load ground-truth origin mappings into dict 
    with open(reads, 'U') as reads_fp:
        for label, seq in parse_fasta(reads_fp):
            expected_reads += 1
            name = label.split()[0]
            length = len(seq)
            if name not in expected_mapping:
                expected_mapping[name] = [length]
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
            sw_score = line.split()[11].split(':')[2]
            # segemehl splits by tabs (label can have spaces)
            contig = line.split('\t')[2]
            if contig == "*":
                continue
            orig_begin = int(line.split('\t')[3])
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
                s1_mapping[name] = [contig, orig_begin, cigar, align_len, num_indel, sw_score]
            else:
                raise ValueError("%s: %s seen twice" % (software_1, name))
    s2_mapping = {}

    print "Done reading s1 file."

    # software 2
    with open(s2_mapping_f, 'U') as sam_fp:
        for line in sam_fp:
            if line.startswith('@'):
                continue
            name = line.split()[0]
            sw_score = line.split()[11].split(':')[2]
            contig = line.split('\t')[2]
            if contig == "*":
                continue
            orig_begin = int(line.split('\t')[3])
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
            if name not in s2_mapping:
                s2_mapping[name] = [contig, orig_begin, cigar, align_len, num_indel, sw_score]
            else:
                raise ValueError("%s: %s seen twice" % (software_2, name))
            
    print "Done reading s2 file."

    # compare ground-truth to both software
    print "Number %s alignments: %s" % (software_1, len(s1_mapping))
    print "Number %s alignments: %s" % (software_2, len(s2_mapping))
    print "read\texpected\t%s\t%s" % (software_1, software_2)
    
    same_alignments = 0
    sw_score_1_higher = 0
    sw_score_2_higher = 0
    cigar_1_longer = 0
    cigar_2_longer = 0
    common_reads = 0
    sw_score_equal = 0
    total_indel_same = 0
    total_read_length = 0
    total_query_coverage = 0.0
    different_scores = 0
    reads_by_1 = 0
    total_read_length_1 = 0
    total_q_cov_1 = 0.0
    total_indel = 0

    for name in expected_mapping:
        if ((name in s1_mapping) and (name in s2_mapping)):
            common_reads += 1
            contig_1 = s1_mapping[name][0]
            orig_begin_1 = s1_mapping[name][1]
            cigar_1 = s1_mapping[name][2]
            cigar_len_1 = s1_mapping[name][3]
            sw_score_1 = s1_mapping[name][5]
            num_indel_1 = s1_mapping[name][4]
            contig_2 = s2_mapping[name][0]
            orig_begin_2 = s2_mapping[name][1]
            cigar_2 = s2_mapping[name][2]
            cigar_len_2 = s2_mapping[name][3]
            sw_score_2 = s2_mapping[name][5]
            num_indel_2 = s2_mapping[name][4]

            # compare cigar lengths
            if (cigar_len_1 > cigar_len_2):
                cigar_1_longer += 1
            if (cigar_len_1 < cigar_len_2):
                cigar_2_longer += 1

            # compare alignments
            if ((contig_1 == contig_2) and (abs(orig_begin_1 - orig_begin_2) <= offset)):
                same_alignments +=1
                if sw_score_1 != sw_score_2:
                    different_scores +=1
                #if num_indel_1 != num_indel_2:
                #    print "num_indel_1 = %s\tnum_indel_2 = %s" % (num_indel_1, num_indel_2)
                #    print "AS_1 = %s\tAS_2 = %s" % (sw_score_1, sw_score_2)
                #    print "cigar_1 = %s\tcigar_2 = %s" % (cigar_1, cigar_2)
                #print "read_len = %s\tcigar_len_1 = %s, query_cov = %s\tcigar_len_2 = %s, query_cov = %s" % (expected_mapping[name][0], cigar_len_1, float(float(cigar_len_1)/float(expected_mapping[name][0])), cigar_len_2, float(float(cigar_len_2)/float(expected_mapping[name][0]))) 
                total_indel_same += num_indel_1
                total_read_length += expected_mapping[name][0]
                total_query_coverage += float(float(cigar_len_1)/float(expected_mapping[name][0]))
            else:
                if (sw_score_1 > sw_score_2):
                    sw_score_1_higher += 1
                    #print "%s\tcontig=%s,pos=%s,cigar=%s,len=%s,sw=%s\tcontig=%s,pos=%s,cigar=%s,len=%s,sw=%s" % (name, contig_1, orig_begin_1, cigar_1, cigar_len_1, sw_score_1, contig_2, orig_begin_2, cigar_2, cigar_len_2, sw_score_2)
                elif (sw_score_1 < sw_score_2):
                    sw_score_2_higher += 1
                elif (sw_score_1 == sw_score_2):
                    sw_score_equal += 1

            #print "%s\tcontig=%s,pos=%s,cigar=%s,len=%s\tcontig=%s,pos=%s,cigar=%s,len=%s" % (name, contig_1, orig_begin_1, cigar_1, cigar_len_1, contig_2, orig_begin_2, cigar_2, cigar_len_2)

        elif ((name in s1_mapping) and (name not in s2_mapping)):
            reads_by_1 += 1
            total_indel += s1_mapping[name][4]
            total_read_length_1 += expected_mapping[name][0]
            total_q_cov_1 += float(float(s1_mapping[name][3])/float(expected_mapping[name][0]))

    print "Number of common reads mapped = %s" % common_reads
    print "Number of same alignments = %s" % same_alignments
    print "Number of different scores = %s" % different_scores
    print "Average numder of indels per alignment: %s" % float(total_indel_same/same_alignments)
    print "Average read length: %s" % int(total_read_length/same_alignments)
    print "Average query coverage: %s" % float(float(total_query_coverage)/float(same_alignments))
    print "Number of longer CIGAR strings by %s = %s" % (software_1, cigar_1_longer)
    print "Number of longer CIGAR strings by %s = %s" % (software_2, cigar_2_longer)
    print "For different alignments, number of alignments with higher SW score for %s = %s" % (software_1, sw_score_1_higher)
    print "For different alignments, number of alignments with higher SW score for %s = %s" % (software_2, sw_score_2_higher)
    print "For different alignments, number of alignments with equal SW score = %s" % sw_score_equal

    print ""
    print "Total reads found only by %s = %s" % (software_1, reads_by_1)
    print "Average number of indels per alignment: %s" % float(total_indel/reads_by_1)
    print "Average read length: %s" % int(total_read_length_1/reads_by_1)
    print "Average query coverage: %s" % float(float(total_q_cov_1)/float(reads_by_1))
