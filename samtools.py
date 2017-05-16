#!/usr/bin/env python

"""
The hilic samtools utilities
"""

import time
import sys
import os
import utils
import subprocess


def run_samtools_concatenation(bam_file_list, combined_bam_filename):
    """
    Concatenate a list of bam file into a single bam file without disturbing the orignial order
    :param bam_file_list: a list of bam files
    :param combined_bam_filename: output combined bam file name
    """
    args = ["samtools", "cat", "-o", combined_bam_filename]
    for filename in bam_file_list:
        args.append(filename)
    subprocess.call(args)


def run_samtools_sort(input_bam_filename, output_bam_prefix, threads):
    """
    Sort bam file by position
    :param input_bam_filename: input bam file name.
    :param output_bam_prefix: output bam file prefix
    :param threads: number of threads for sorting
    """
    args = ["samtools", "sort", "-@", str(threads), input_bam_filename, "-o", output_bam_prefix]
    subprocess.call(args)


def run_samtools_index(sorted_bam_filename):
    """
    Index the sorted bam file
    :param sorted_bam_filename: input bam file to index.
    """
    args = ["samtools", "index", sorted_bam_filename]
    subprocess.call(args)
