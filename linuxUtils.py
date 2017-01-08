#!/usr/bin/env python

"""
The Linux wrapper to run command line programs in hilic.py
"""

import subprocess
import sys
import os


def run_bwa(fastq_filename, genome_fasta, sam_filename):
    DEVNULL = open(os.devnull, 'w')
    f = open(sam_filename, 'w')
    args = ["bwa", "mem", "-5", genome_fasta, fastq_filename]
    process = subprocess.Popen(args, stdout=f, stderr=DEVNULL)
    return (process, fastq_filename)



