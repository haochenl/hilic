#!/usr/bin/env python

"""
The Linux wrapper to run command line programs in hilic.py
"""

import os
import subprocess


def run_bwa(fastq_filename, genome_fasta, sam_filename):
    DEVNULL = open(os.devnull, 'w')
    f = open(sam_filename, 'w')
    args = ["bwa", "mem", "-5", genome_fasta, fastq_filename]
    process = subprocess.Popen(args, stdout=f, stderr=DEVNULL)
    return (process, fastq_filename)


def run_samtools(sam_filename, bam_filename):
    DEVNULL = open(os.devnull, 'w')
    f = open(bam_filename, 'w')
    p1_args = ["samtools", "view", "-S", "-F", "2048", "-h", sam_filename]
    p1 = subprocess.Popen(p1_args, stdout=subprocess.PIPE, stderr=DEVNULL)
    p2_args = ["samtools", "view", "-bS", "-"]
    p2 = subprocess.Popen(p2_args, stdin=p1.stdout, stdout=f, stderr=DEVNULL)
    p1.stdout.close()
    return (p2, sam_filename)
