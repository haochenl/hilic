#!/usr/bin/env python

"""
The Linux wrapper to run command line programs in hilic.py
"""

import subprocess
import sys


def run_bwa(fastq_filename, genome_fasta, sam_filename):
    args = ["bwa", "mem", "-5", genome_fasta, fastq_filename, ">", sam_filename]
    code = subprocess.call(args)
    if code == 0:
        print >> sys.stderr, "File '{}' processed.".format(fastq_filename)
    else:
        print >> sys.stderr, "Error when processing file '{}'".format(fastq_filename)
