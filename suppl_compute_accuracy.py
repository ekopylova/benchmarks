#!/opt/python-2.7.3/bin/python

"""Compute the accuracy of a tool's alignments vs. ground-truth alignments
   usage: python suppl_compute_accuracy.py ground_truth.blast tool [sam,blast] tool_alignment.[sam,blast] total_reads
"""

import sys
import click

def collect_ground_truth(ground_truth_alns_fp):
    """Parses a file of ground truth alignments into a dictionary
       (only BLAST)

       Parameters:
       -----------
       ground_truth_alns_fp : string
          filepath of ground-truth BLAST alignments

       Returns:
       --------
       expected_alns : dict
          dictionary of ground-truth alignments
    """
    expected_alns = {}
    # collect ground-truth alignments
    with open(ground_truth_alns_fp, 'U') as ground_truth_alns:
        for line in ground_truth_alns:
            line = line.strip().split()
            read_id = line[0]
            if read_id in expected_alns:
                expected_alns[read_id].append(line[1:])
            else:
                expected_alns[read_id] = [line[1:]]

    # sort alignments for each read by highest bitscore
    for read in expected_alns:
        expected_alns[read].sort(key=lambda x: float(x[10]), reverse=True)

    return expected_alns


def collect_observed_alignments(observed_alns_fp, file_format="sam"):
    """Parses a file of observed alignments into a dictionary
       (BLAST or SAM)

       Parameters:
       -----------
       observed_alns_fp : string
          filepath to observed alignments
       file_format : string
          file format of alignments (blast or sam)

       Returns:
       --------
       observed_alns : dict
          dictionary of observed alignments
    """
    observed_alns = {}

    with open(observed_alns_fp, 'U') as observed_alns_f:
        if file_format == "sam":
            for line in observed_alns_f:
                if line.startswith('@'):
                    continue
                line = line.strip().split('\t')
                contig = line[2].split()
                if contig:
                    contig = contig[0]
                # no alignment found for this read
                if contig == "*":
                    continue

                read_id = line[0].split()
                if read_id:
                    read_id = read_id[0]
                if read_id not in observed_alns:
                    observed_alns[read_id] = line[1:]
                else:
                    raise ValueError("Only 1 alignment per read: %s" % read_id)

        elif file_format == "blast":
            for line in observed_alns_f:
                line = line.strip().split('\t')
                contig = line[1]
                read_id = line[0]
                if read_id not in observed_alns:
                    observed_alns[read_id] = line[1:]
                else:
                    raise ValueError("Only 1 alignment per read: %s" % read_id)
        else:
            raise ValueError("%s file format not supported" % file_format)

    return observed_alns


def compute_accuracy(expected_alns, observed_alns, file_format="sam", offset=0):
    """For each observed alignment, compute the accuracy score based on the list of
       expected alignments. The accuracy score is between [0,1] and is weighted based
       on the quality of the alignment given by the bit score.

       For example, given the following 6 expected alignments for a read:
          1. --------- bitscore = 400 (score 4/4 = 1)
          2. --------- bitscore = 400 (score 4/4 = 1)
          3. --------- bitscore = 350 (score 3/4 = 0.75)
          4. --------- bitscore = 350 (score 3/4 = 0.75)
          5. --------- bitscore = 250 (score 2/4 = 0.5)
          6. --------- bitscore = 150 (score 1/4 = 0.25)
       If the observed alignment for the same read matches by contig and position
       to expected alignment # 1 or # 2, it will have an accuracy score of 1.

       Parameters:
       -----------
       expected_alns : dictionary
          dictionary of expected alignments, keys are read ids and values are
          a list of of BLAST alignments
       observed_alns : dictionary
          dictionary of observed alignments, keys are read ids and values are
          a list of SAM or BLAST alignments
       file_format : string, optional
          file format of observed alignments (SAM or BLAST)
       offset : integer, optional
          the maximum difference between the expected alignment position and
          the observed

       Returns:
       --------
       total_accuracy_score : float
          the total accuracy score between [0,1] for all alignments in observed
          alignments
    """
    total_accuracy_score = 0.0
    all_accuracy_scores = []
    false_negative = set()
    false_positive = set()
    true_positive = set()

    if file_format == "sam":
        obs_contig_index = 1
        obs_pos_index = 2
    elif file_format == "blast":
        obs_contig_index = 0
        obs_pos_index = 7
    else:
        raise ValueError("%s file format not supported" % file_format)

    for read_id in observed_alns:
        if read_id not in expected_alns:
            #print "WARNING: alignment %s in observed but not expected" % read_id
            continue
        obs_contig = observed_alns[read_id][obs_contig_index]
        obs_pos = int(observed_alns[read_id][obs_pos_index])
        if file_format == "blast":
            # read mapped as reverse-complement
            if (int(observed_alns[read_id][7]) - int(observed_alns[read_id][8]) > 0):
                obs_strand = 16
            else:
                obs_strand = 0
        # compute the total number of unique bitscores in list of expected alignments
        unique_bitscores = set()
        for exp_aln in expected_alns[read_id]:
            unique_bitscores.add(exp_aln[10])
        num_unique_bitscores = float(len(unique_bitscores))
        # find the expected alignment for this read
        index = 0
        weight = num_unique_bitscores
        for exp_aln in expected_alns[read_id]:
            exp_contig = exp_aln[0]
            if ((file_format == "blast") and (obs_strand == 16)):
                exp_pos = int(exp_aln[8])
            else:
                exp_pos = int(exp_aln[7])
            # alignment found, compute accuracy score and go to next read
            if ((obs_contig == exp_contig) and (abs(obs_pos - exp_pos) <= offset)):
                accuracy_score = float(weight/num_unique_bitscores)
                all_accuracy_scores.append(accuracy_score)
                break
            if index < len(expected_alns[read_id]) - 1:
                # decrement the weight
                if float(expected_alns[read_id][index+1][10]) < float(exp_aln[10]):
                    weight -= 1
            index += 1

    total_accuracy_score = float(sum(all_accuracy_scores)/len(observed_alns))
        
    return total_accuracy_score*100.0


def compute_precision(expected_alns, observed_alns):
    """
    """
    expected_alns_keys = set(expected_alns.keys())
    observed_alns_keys = set(observed_alns.keys())

    # compute true positive, false positive and false negative read counts
    tp = len(observed_alns_keys & expected_alns_keys)
    fp = len(observed_alns_keys - expected_alns_keys)
    fn = len(expected_alns_keys - observed_alns_keys)

    # compute precision, recall and F-measure for read counts
    p = tp / float(tp + fp)
    r = tp / float(tp + fn)
    f = float(2 * p * r) / float(p + r)

    return tp, fp, fn, p, r, f


@click.command()
@click.argument('expected_alns_fp', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=True,
                                file_okay=True))
@click.argument('observed_alns_fp', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=True,
                                file_okay=True))
@click.option('--tool', type=str, required=True,
              help='tool which generated observed_alns_fp')
@click.option('--offset', type=int, required=False,
              help='maximum absolute difference between expected and observed origin position')
@click.option('--observed_aln_format', type=str, default='sam', show_default=True,
              required=False, help='file format of observed_alns_fp')
@click.option('--total_reads', type=int, required=True,
              help='total number of sequences used in alignment for generating observed_alns_fp')
def _main(expected_alns_fp, observed_alns_fp, tool, offset, observed_aln_format, total_reads):
    """
    """
    allowed_alignment_types = ["sam", "blast"]
    if observed_aln_format not in allowed_alignment_types:
        raise ValueError("%s is not supported" % observed_aln_format)

    expected_alns = collect_ground_truth(ground_truth_alns_fp=expected_alns_fp)
    observed_alns = collect_observed_alignments(observed_alns_fp=observed_alns_fp, file_format=observed_aln_format)
    accuracy = compute_accuracy(expected_alns=expected_alns,
                                observed_alns=observed_alns,
                                file_format=observed_aln_format,
                                offset=offset)

    tp, fp, fn, p, r, f = compute_precision(expected_alns, observed_alns)

    total_reads_mapped = float(len(observed_alns))/float(total_reads)

    sys.stdout.write("%s\t%.2f\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.2f\t" %
        (len(observed_alns), total_reads_mapped, tp, fp, fn, p, r, f, accuracy))


if __name__ == "__main__":
    _main()
