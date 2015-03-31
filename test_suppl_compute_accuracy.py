#!/usr/bin/env python
"""
Unit tests for suppl_compute_accuracy.py
========================================
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, Evguenia Kopylova
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

from skbio.util import remove_files
from tempfile import mkstemp
from os import close

from suppl_compute_accuracy import (collect_ground_truth,
                                    collect_observed_alignments,
                                    compute_accuracy,
                                    compute_precision)


# Test class and cases
class SupplComputeAccuracyTests(TestCase):
    """ Tests for rsuppl_compute_accuracy.py functionality """

    def setUp(self):
        """ 
        """
        # write expected_alignments_1 to temporary file
        f, self.exp_alns_1_fp = mkstemp(prefix='exp_alns_1_',
                                        suffix='.txt')
        close(f)
        with open(self.exp_alns_1_fp, 'w') as tmp:
            tmp.write(expected_alignments_1)

        # write observed_sam_alignments_1 to temporary file
        f, self.obs_sam_alns_1_fp = mkstemp(prefix='obs_alns_1_',
                                            suffix='.sam')
        close(f)
        with open(self.obs_sam_alns_1_fp, 'w') as tmp:
            tmp.write(observed_sam_alignments_1)

        # write observed_blast_alignments_1 to temporary file
        f, self.obs_blast_alns_1_fp = mkstemp(prefix='obs_alns_1_',
                                              suffix='.blast')
        close(f)
        with open(self.obs_blast_alns_1_fp, 'w') as tmp:
            tmp.write(observed_blast_alignments_1)

        self.files_to_remove = [self.exp_alns_1_fp,
                                self.obs_sam_alns_1_fp,
                                self.obs_blast_alns_1_fp]

    def tearDown(self):
        remove_files(self.files_to_remove)

    def test_collect_ground_truth(self):
        """
        """
        exp_alns = {
            'seq.000000828': [
            ['ref1', '96.34', '82', '3', '0', '69', '150', '2428038', '2428119', '9.5e-26', '111.8'],
            ['ref1', '98.61', '72', '1', '0', '1', '72', '2426641', '2426712', '1.6e-23', '104.5']],
            'seq.000001026': [
            ['ref1', '98.00', '150', '0', '3', '1', '150', '524183', '524329', '3.5e-69', '256.1'],
            ['ref1', '97.33', '150', '1', '3', '1', '150', '568163', '568309', '8.1e-68', '251.6'],
            ['ref1', '97.33', '150', '1', '3', '150', '1', '2292122', '2292268', '8.1e-68', '251.6'],
            ['ref1', '97.33', '150', '1', '3', '150', '1', '2174250', '2174396', '8.1e-68', '251.6'],
            ['ref1', '96.67', '150', '2', '3', '150', '1', '1978335', '1978481', '1.8e-66', '247.1']],
            'seq.000001004': [
            ['ref1', '99.33', '150', '1', '0', '150', '1', '2290660', '2290809', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '150', '1', '1976873', '1977022', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '1', '150', '525642', '525791', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '1', '150', '569622', '569771', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '150', '1', '2172788', '2172937', '3.9e-74', '272.6']]}

        given_alns = collect_ground_truth(self.exp_alns_1_fp)

        self.assertEqual(len(exp_alns), len(given_alns))

        for aln in given_alns:
            self.assertTrue(aln in exp_alns)
            # list of lists must be sorted by decreasing bit score
            self.assertEqual(exp_alns[aln], given_alns[aln])

    def test_collect_observed_sam_alignments(self):
        """
        """
        exp_sam_alns = {
            'seq.000000828':
            ['0', 'ref1', '2428038', '255', '68S82M', '*', '0', '0', 'GATTAATTTATAAAAATCAAATCAGGCATTAAATAAAATAGCCCATAAATACAAAGTGTTATCACCTTCTATTTAGGGGCTATTAGTTCTATTCGTCATTCTTTTTACAGATCATTCTATCTAATTAATTTGTGTACAATTTTGATAACT', '*', 'AS:i:149', 'NM:i:3'],
            'seq.000001026':
            ['0', 'ref1', '524183', '255', '110M3I37M', '*', '0', '0', 'CCCGGAGGAAGAGAAAGAAAATTCGATTCCCTTAGTAGCGGCGAGCGAAATGGGAAGAGCCCAAACCAACAAGCTTGCTTGTTGGGGTTGTAGGACACTCTATACGGAGTCTCTACAAAGGACGACATTAGACGAATCATCTGGAAAGAT', '*', 'AS:i:285', 'NM:i:3'],
            'seq.000001004':
            ['0', 'ref1', '525642', '255', '150M', '*', '0', '0', 'AAGATGAGAATTCTAAGGTGAGCGAGCGAACTCTCGTTAAGGAACTCGGCAAAATGACCCCGTAACTTCGGGAGAAGGGGTGCTCTTTAGGGTTAACGCCCAGAAGAGCCGCAGTGAATAGGCCCAAGCGACTGTTTATCAAACACACAG', '*', 'AS:i:295', 'NM:i:1']}

        obs_alns = collect_observed_alignments(self.obs_sam_alns_1_fp)

        self.assertEqual(len(exp_sam_alns), len(obs_alns))

        for aln in obs_alns:
            self.assertTrue(aln in exp_sam_alns)
            self.assertEqual(exp_sam_alns[aln], obs_alns[aln])

    def test_collect_observed_blast_alignments(self):
        """
        """
        exp_blast_alns = {
        'seq.000000828': ['ref1', '96.34', '82', '3', '0', '69', '150', '2428038', '2428119', '1e-32', '134'],
        'seq.000001026': ['ref1', '98.00', '150', '0', '1', '1', '150', '524183', '524329', '4e-69', '255'],
        'seq.000001004': ['ref1', '99.33', '150', '1', '0', '1', '150', '525642', '525791', '2e-72', '266']
        }

        obs_alns = collect_observed_alignments(observed_alns_fp=self.obs_blast_alns_1_fp,
                                               file_format="blast")

        self.assertEqual(len(exp_blast_alns), len(obs_alns))

        for aln in obs_alns:
            self.assertTrue(aln in exp_blast_alns)
            self.assertEqual(exp_blast_alns[aln], obs_alns[aln])      

    def test_compute_accuracy_100_sam(self):
        """ Given a set of expected alignments, each observed alignment
            should achieve an accuracy score of 100.0%, thus a total of
            100.0% accuracy score for all observed alignments
            (observed alignments in SAM format)
        """
        exp_alns = {
            'seq.000000828': [
            ['ref1', '96.34', '82', '3', '0', '69', '150', '2428038', '2428119', '9.5e-26', '111.8'],
            ['ref1', '98.61', '72', '1', '0', '1', '72', '2426641', '2426712', '1.6e-23', '104.5']],
            'seq.000001026': [
            ['ref1', '98.00', '150', '0', '3', '1', '150', '524183', '524329', '3.5e-69', '256.1'],
            ['ref1', '97.33', '150', '1', '3', '1', '150', '568163', '568309', '8.1e-68', '251.6'],
            ['ref1', '97.33', '150', '1', '3', '150', '1', '2292122', '2292268', '8.1e-68', '251.6'],
            ['ref1', '97.33', '150', '1', '3', '150', '1', '2174250', '2174396', '8.1e-68', '251.6'],
            ['ref1', '96.67', '150', '2', '3', '150', '1', '1978335', '1978481', '1.8e-66', '247.1']],
            'seq.000001004': [
            ['ref1', '99.33', '150', '1', '0', '150', '1', '2290660', '2290809', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '150', '1', '1976873', '1977022', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '1', '150', '525642', '525791', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '1', '150', '569622', '569771', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '150', '1', '2172788', '2172937', '3.9e-74', '272.6']]}

        obs_alns = {
            'seq.000000828':
            ['0', 'ref1', '2428038', '255', '68S82M', '*', '0', '0', 'GATTAATTTATAAAAATCAAATCAGGCATTAAATAAAATAGCCCATAAATACAAAGTGTTATCACCTTCTATTTAGGGGCTATTAGTTCTATTCGTCATTCTTTTTACAGATCATTCTATCTAATTAATTTGTGTACAATTTTGATAACT', '*', 'AS:i:149', 'NM:i:3'],
            'seq.000001026':
            ['0', 'ref1', '524183', '255', '110M3I37M', '*', '0', '0', 'CCCGGAGGAAGAGAAAGAAAATTCGATTCCCTTAGTAGCGGCGAGCGAAATGGGAAGAGCCCAAACCAACAAGCTTGCTTGTTGGGGTTGTAGGACACTCTATACGGAGTCTCTACAAAGGACGACATTAGACGAATCATCTGGAAAGAT', '*', 'AS:i:285', 'NM:i:3'],
            'seq.000001004':
            ['0', 'ref1', '525642', '255', '150M', '*', '0', '0', 'AAGATGAGAATTCTAAGGTGAGCGAGCGAACTCTCGTTAAGGAACTCGGCAAAATGACCCCGTAACTTCGGGAGAAGGGGTGCTCTTTAGGGTTAACGCCCAGAAGAGCCGCAGTGAATAGGCCCAAGCGACTGTTTATCAAACACACAG', '*', 'AS:i:295', 'NM:i:1']}

        accuracy = compute_accuracy(exp_alns, obs_alns)

        self.assertEqual(accuracy, 100.0)

    def test_compute_accuracy_100_blast(self):
        """ Given a set of expected alignments, each observed alignment
            should achieve an accuracy score of 100.0%, thus a total of
            100.0% accuracy score for all observed alignments
            (observed alignments in BLAST format)
        """
        exp_alns = {
            'seq.000000828': [
            ['ref1', '96.34', '82', '3', '0', '69', '150', '2428038', '2428119', '9.5e-26', '111.8'],
            ['ref1', '98.61', '72', '1', '0', '1', '72', '2426641', '2426712', '1.6e-23', '104.5']],
            'seq.000001026': [
            ['ref1', '98.00', '150', '0', '3', '1', '150', '524183', '524329', '3.5e-69', '256.1'],
            ['ref1', '97.33', '150', '1', '3', '1', '150', '568163', '568309', '8.1e-68', '251.6'],
            ['ref1', '97.33', '150', '1', '3', '150', '1', '2292122', '2292268', '8.1e-68', '251.6'],
            ['ref1', '97.33', '150', '1', '3', '150', '1', '2174250', '2174396', '8.1e-68', '251.6'],
            ['ref1', '96.67', '150', '2', '3', '150', '1', '1978335', '1978481', '1.8e-66', '247.1']],
            'seq.000001004': [
            ['ref1', '99.33', '150', '1', '0', '150', '1', '2290660', '2290809', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '150', '1', '1976873', '1977022', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '1', '150', '525642', '525791', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '1', '150', '569622', '569771', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '150', '1', '2172788', '2172937', '3.9e-74', '272.6']]}

        obs_alns = {
        'seq.000000828': ['ref1', '96.34', '82', '3', '0', '69', '150', '2428038', '2428119', '1e-32', '134'],
        'seq.000001026': ['ref1', '98.00', '150', '0', '1', '1', '150', '524183', '524329', '4e-69', '255'],
        'seq.000001004': ['ref1', '99.33', '150', '1', '0', '1', '150', '525642', '525791', '2e-72', '266']
        }

        accuracy = compute_accuracy(exp_alns, obs_alns, file_format="blast")

        self.assertEqual(accuracy, 100.0)

    def test_compute_accuracy_72_sam(self):
        """ Given a set of expected alignments, 2/3 observed alignments
            are not optimal with seq.000000828 & seq.000001026 being second best,
            and seq.000001004 being best.
            The final accuracy score is (1/2 + 2/3 + 1)/3 = 0.722
            (observed alignments in SAM format)
        """
        exp_alns = {
            'seq.000000828': [
            ['ref1', '96.34', '82', '3', '0', '69', '150', '2428038', '2428119', '9.5e-26', '111.8'],
            ['ref1', '98.61', '72', '1', '0', '1', '72', '2426641', '2426712', '1.6e-23', '104.5']],
            'seq.000001026': [
            ['ref1', '98.00', '150', '0', '3', '1', '150', '524183', '524329', '3.5e-69', '256.1'],
            ['ref1', '97.33', '150', '1', '3', '1', '150', '568163', '568309', '8.1e-68', '251.6'],
            ['ref1', '97.33', '150', '1', '3', '150', '1', '2292122', '2292268', '8.1e-68', '251.6'],
            ['ref1', '97.33', '150', '1', '3', '150', '1', '2174250', '2174396', '8.1e-68', '251.6'],
            ['ref1', '96.67', '150', '2', '3', '150', '1', '1978335', '1978481', '1.8e-66', '247.1']],
            'seq.000001004': [
            ['ref1', '99.33', '150', '1', '0', '150', '1', '2290660', '2290809', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '150', '1', '1976873', '1977022', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '1', '150', '525642', '525791', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '1', '150', '569622', '569771', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '150', '1', '2172788', '2172937', '3.9e-74', '272.6']]}

        obs_alns = {
            'seq.000000828':
            ['0', 'ref1', '2426641', '255', '72M78S', '*', '0', '0', 'GATTAATTTATAAAAATCAAATCAGGCATTAAATAAAATAGCCCATAAATACAAAGTGTTATCACCTTCTATTTAGGGGCTATTAGTTCTATTCGTCATTCTTTTTACAGATCATTCTATCTAATTAATTTGTGTACAATTTTGATAACT', '*', 'AS:i:139', 'NM:i:1'],
            'seq.000001026':
            ['0', 'ref1', '568163', '255', '110M3I37M', '*', '0', '0', 'CCCGGAGGAAGAGAAAGAAAATTCGATTCCCTTAGTAGCGGCGAGCGAAATGGGAAGAGCCCAAACCAACAAGCTTGCTTGTTGGGGTTGTAGGACACTCTATACGGAGTCTCTACAAAGGACGACATTAGACGAATCATCTGGAAAGAT', '*', 'AS:i:280', 'NM:i:4'],
            'seq.000001004':
            ['0', 'ref1', '525642', '255', '150M', '*', '0', '0', 'AAGATGAGAATTCTAAGGTGAGCGAGCGAACTCTCGTTAAGGAACTCGGCAAAATGACCCCGTAACTTCGGGAGAAGGGGTGCTCTTTAGGGTTAACGCCCAGAAGAGCCGCAGTGAATAGGCCCAAGCGACTGTTTATCAAACACACAG', '*', 'AS:i:295', 'NM:i:1']}

        accuracy = compute_accuracy(exp_alns, obs_alns)

        self.assertEqual(float("%0.1f" % accuracy), 72.2)

    def test_compute_accuracy_72_blast(self):
        """ Given a set of expected alignments, 2/3 observed alignments
            are not optimal with seq.000000828 & seq.000001026 being second best,
            and seq.000001004 being best.
            The final accuracy score is (1/2 + 2/3 + 1)/3 = 0.722
            (observed alignments in BLAST format)
        """
        exp_alns = {
            'seq.000000828': [
            ['ref1', '96.34', '82', '3', '0', '69', '150', '2428038', '2428119', '9.5e-26', '111.8'],
            ['ref1', '98.61', '72', '1', '0', '1', '72', '2426641', '2426712', '1.6e-23', '104.5']],
            'seq.000001026': [
            ['ref1', '98.00', '150', '0', '3', '1', '150', '524183', '524329', '3.5e-69', '256.1'],
            ['ref1', '97.33', '150', '1', '3', '1', '150', '568163', '568309', '8.1e-68', '251.6'],
            ['ref1', '97.33', '150', '1', '3', '150', '1', '2292122', '2292268', '8.1e-68', '251.6'],
            ['ref1', '97.33', '150', '1', '3', '150', '1', '2174250', '2174396', '8.1e-68', '251.6'],
            ['ref1', '96.67', '150', '2', '3', '150', '1', '1978335', '1978481', '1.8e-66', '247.1']],
            'seq.000001004': [
            ['ref1', '99.33', '150', '1', '0', '150', '1', '2290660', '2290809', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '150', '1', '1976873', '1977022', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '1', '150', '525642', '525791', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '1', '150', '569622', '569771', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '150', '1', '2172788', '2172937', '3.9e-74', '272.6']]}

        obs_alns = {
            'seq.000000828':
            ['ref1', '98.61', '72', '1', '0', '1', '72', '2426641', '2426712', '5e-30', '125'],
            'seq.000001026':
            ['ref1', '97.33', '150', '1', '1', '1', '150', '568163', '568309', '5e-68', '251'],
            'seq.000001004':
            ['ref1', '99.33', '150', '1', '0',  '1', '150', '525642', '525791', '2e-72', '266']}

        accuracy = compute_accuracy(exp_alns, obs_alns, file_format="blast")

        self.assertEqual(float("%0.1f" % accuracy), 72.2)

    def test_compute_precision_100(self):
        """Test functionality of compute_precision() method,
           expected to return 100% precision, recall and F-measure
        """
        exp_alns = {
            'seq.000000828': [
            ['ref1', '96.34', '82', '3', '0', '69', '150', '2428038', '2428119', '9.5e-26', '111.8'],
            ['ref1', '98.61', '72', '1', '0', '1', '72', '2426641', '2426712', '1.6e-23', '104.5']],
            'seq.000001026': [
            ['ref1', '98.00', '150', '0', '3', '1', '150', '524183', '524329', '3.5e-69', '256.1'],
            ['ref1', '97.33', '150', '1', '3', '1', '150', '568163', '568309', '8.1e-68', '251.6'],
            ['ref1', '97.33', '150', '1', '3', '150', '1', '2292122', '2292268', '8.1e-68', '251.6'],
            ['ref1', '97.33', '150', '1', '3', '150', '1', '2174250', '2174396', '8.1e-68', '251.6'],
            ['ref1', '96.67', '150', '2', '3', '150', '1', '1978335', '1978481', '1.8e-66', '247.1']],
            'seq.000001004': [
            ['ref1', '99.33', '150', '1', '0', '150', '1', '2290660', '2290809', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '150', '1', '1976873', '1977022', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '1', '150', '525642', '525791', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '1', '150', '569622', '569771', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '150', '1', '2172788', '2172937', '3.9e-74', '272.6']]}

        obs_alns = {
            'seq.000000828':
            ['0', 'ref1', '2426641', '255', '72M78S', '*', '0', '0', 'GATTAATTTATAAAAATCAAATCAGGCATTAAATAAAATAGCCCATAAATACAAAGTGTTATCACCTTCTATTTAGGGGCTATTAGTTCTATTCGTCATTCTTTTTACAGATCATTCTATCTAATTAATTTGTGTACAATTTTGATAACT', '*', 'AS:i:139', 'NM:i:1'],
            'seq.000001026':
            ['0', 'ref1', '568163', '255', '110M3I37M', '*', '0', '0', 'CCCGGAGGAAGAGAAAGAAAATTCGATTCCCTTAGTAGCGGCGAGCGAAATGGGAAGAGCCCAAACCAACAAGCTTGCTTGTTGGGGTTGTAGGACACTCTATACGGAGTCTCTACAAAGGACGACATTAGACGAATCATCTGGAAAGAT', '*', 'AS:i:280', 'NM:i:4'],
            'seq.000001004':
            ['0', 'ref1', '525642', '255', '150M', '*', '0', '0', 'AAGATGAGAATTCTAAGGTGAGCGAGCGAACTCTCGTTAAGGAACTCGGCAAAATGACCCCGTAACTTCGGGAGAAGGGGTGCTCTTTAGGGTTAACGCCCAGAAGAGCCGCAGTGAATAGGCCCAAGCGACTGTTTATCAAACACACAG', '*', 'AS:i:295', 'NM:i:1']}

        tp, fp, fn, p, r, f = compute_precision(exp_alns, obs_alns)

        self.assertEqual(tp, 3)
        self.assertEqual(fn, 0)
        self.assertEqual(fp, 0)
        self.assertEqual(p, 1.0)
        self.assertEqual(r, 1.0)
        self.assertEqual(f, 1.0)

    def test_compute_precision_67(self):
        """Test functionality of compute_precision() method,
           expected to return 67% precision, 40% recall and 50% F-measure
        """
        exp_alns = {
            'seq.000000000': [
            ['ref1', '95.60', '273', '12', '0', '274', '2', '1409508', '1409780', '9.9e-111', '395.0']],
            'seq.000000003': [
            ['ref1', '98.18', '275', '4', '1', '275', '1', '144640', '144913', '1.8e-116', '414.1'],
            ['ref1', '73.33', '150', '28', '12', '126', '269', '2859998', '2860141', '1e-10', '62.7']],
            'seq.000000828': [
            ['ref1', '96.34', '82', '3', '0', '69', '150', '2428038', '2428119', '9.5e-26', '111.8'],
            ['ref1', '98.61', '72', '1', '0', '1', '72', '2426641', '2426712', '1.6e-23', '104.5']],
            'seq.000001026': [
            ['ref1', '98.00', '150', '0', '3', '1', '150', '524183', '524329', '3.5e-69', '256.1'],
            ['ref1', '97.33', '150', '1', '3', '1', '150', '568163', '568309', '8.1e-68', '251.6'],
            ['ref1', '97.33', '150', '1', '3', '150', '1', '2292122', '2292268', '8.1e-68', '251.6'],
            ['ref1', '97.33', '150', '1', '3', '150', '1', '2174250', '2174396', '8.1e-68', '251.6'],
            ['ref1', '96.67', '150', '2', '3', '150', '1', '1978335', '1978481', '1.8e-66', '247.1']],
            'seq.000001004': [
            ['ref1', '99.33', '150', '1', '0', '150', '1', '2290660', '2290809', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '150', '1', '1976873', '1977022', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '1', '150', '525642', '525791', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '1', '150', '569622', '569771', '3.9e-74', '272.6'],
            ['ref1', '99.33', '150', '1', '0', '150', '1', '2172788', '2172937', '3.9e-74', '272.6']]}

        obs_alns = {
            'seq.000000828':
            ['0', 'ref1', '2426641', '255', '72M78S', '*', '0', '0', 'GATTAATTTATAAAAATCAAATCAGGCATTAAATAAAATAGCCCATAAATACAAAGTGTTATCACCTTCTATTTAGGGGCTATTAGTTCTATTCGTCATTCTTTTTACAGATCATTCTATCTAATTAATTTGTGTACAATTTTGATAACT', '*', 'AS:i:139', 'NM:i:1'],
            'seq.000001026':
            ['0', 'ref1', '568163', '255', '110M3I37M', '*', '0', '0', 'CCCGGAGGAAGAGAAAGAAAATTCGATTCCCTTAGTAGCGGCGAGCGAAATGGGAAGAGCCCAAACCAACAAGCTTGCTTGTTGGGGTTGTAGGACACTCTATACGGAGTCTCTACAAAGGACGACATTAGACGAATCATCTGGAAAGAT', '*', 'AS:i:280', 'NM:i:4'],
            'seq.000012323':
            ['16' 'ref1', '1911165', '255', '35S3M1D9M1I10M1D25M1I2M1D3M1D17M1I5M1D14M2D3M2I9M1I20M2I11M2D15M21S', '*', '0', '0', 'ACCAATATCTTTCTGATATTGATTATTCTTCCCGAATCTATACTTGCTTGGAATTTATATTTTTTATCATCATTTATATAACCTACAATGAAATACCACCCATCGGAGTTTTATGGAAATCTGTCACTTTCATATTTTTATTAATTGTTAAGATTATGCTTAAAATACAAATCGATTCGTTTTTCTTGTGTCTTAATGTATGAGCCTTTA', '*', 'AS:i:69', 'NM:i:47']}

        tp, fp, fn, p, r, f = compute_precision(exp_alns, obs_alns)

        self.assertEqual(tp, 2)
        self.assertEqual(fn, 3)
        self.assertEqual(fp, 1)
        self.assertEqual(float("%0.2f" % p), 0.67)
        self.assertEqual(float("%0.2f" % r), 0.40)
        self.assertEqual(float("%0.2f" % f), 0.50)


expected_alignments_1 = """seq.000000828\tref1\t96.34\t82\t3\t0\t69\t150\t2428038\t2428119\t9.5e-26\t111.8
seq.000000828\tref1\t98.61\t72\t1\t0\t1\t72\t2426641\t2426712\t1.6e-23\t104.5
seq.000001004\tref1\t99.33\t150\t1\t0\t150\t1\t2290660\t2290809\t3.9e-74\t272.6
seq.000001004\tref1\t99.33\t150\t1\t0\t150\t1\t1976873\t1977022\t3.9e-74\t272.6
seq.000001004\tref1\t99.33\t150\t1\t0\t1\t150\t525642\t525791\t3.9e-74\t272.6
seq.000001004\tref1\t99.33\t150\t1\t0\t1\t150\t569622\t569771\t3.9e-74\t272.6
seq.000001004\tref1\t99.33\t150\t1\t0\t150\t1\t2172788\t2172937\t3.9e-74\t272.6
seq.000001026\tref1\t96.67\t150\t2\t3\t150\t1\t1978335\t1978481\t1.8e-66\t247.1
seq.000001026\tref1\t98.00\t150\t0\t3\t1\t150\t524183\t524329\t3.5e-69\t256.1
seq.000001026\tref1\t97.33\t150\t1\t3\t1\t150\t568163\t568309\t8.1e-68\t251.6
seq.000001026\tref1\t97.33\t150\t1\t3\t150\t1\t2292122\t2292268\t8.1e-68\t251.6
seq.000001026\tref1\t97.33\t150\t1\t3\t150\t1\t2174250\t2174396\t8.1e-68\t251.6
"""

observed_sam_alignments_1 = """@HD  VN:1.0  SO:unsorted
@PG ID:program1    VN:1.0  CL:test1
seq.000000828\t0\tref1\t2428038\t255\t68S82M\t*\t0\t0\tGATTAATTTATAAAAATCAAATCAGGCATTAAATAAAATAGCCCATAAATACAAAGTGTTATCACCTTCTATTTAGGGGCTATTAGTTCTATTCGTCATTCTTTTTACAGATCATTCTATCTAATTAATTTGTGTACAATTTTGATAACT\t*\tAS:i:149\tNM:i:3
seq.000001004\t0\tref1\t525642\t255\t150M\t*\t0\t0\tAAGATGAGAATTCTAAGGTGAGCGAGCGAACTCTCGTTAAGGAACTCGGCAAAATGACCCCGTAACTTCGGGAGAAGGGGTGCTCTTTAGGGTTAACGCCCAGAAGAGCCGCAGTGAATAGGCCCAAGCGACTGTTTATCAAACACACAG\t*\tAS:i:295\tNM:i:1
seq.000001026\t0\tref1\t524183\t255\t110M3I37M\t*\t0\t0\tCCCGGAGGAAGAGAAAGAAAATTCGATTCCCTTAGTAGCGGCGAGCGAAATGGGAAGAGCCCAAACCAACAAGCTTGCTTGTTGGGGTTGTAGGACACTCTATACGGAGTCTCTACAAAGGACGACATTAGACGAATCATCTGGAAAGAT\t*\tAS:i:285\tNM:i:3
"""

observed_blast_alignments_1 = """seq.000000828\tref1\t96.34\t82\t3\t0\t69\t150\t2428038\t2428119\t1e-32\t134
seq.000001026\tref1\t98.00\t150\t0\t1\t1\t150\t524183\t524329\t4e-69\t255
seq.000001004\tref1\t99.33\t150\t1\t0\t1\t150\t525642\t525791\t2e-72\t266
"""

if __name__ == '__main__':
    main()