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
