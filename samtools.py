#!/usr/bin/env python

"""
The hilic samtools utilities
"""

import time
import sys
import os
import utils
import subprocess


def compress_from_reader(input_file_reader):
    """
    Compress the sam alignment files from hilic config file reader
    :param input_file_reader: the hilic config file reader
    :return:
    """
    ## check if a sam file is already compressed
    ## useful when HiC and control files are the same
    compressed_file_set = set()
    compression_processes = []
    compression_start_time = time.time()
    compress_sam_list(input_file_reader.hicRead1Files, compressed_file_set, compression_processes)
    compress_sam_list(input_file_reader.hicRead2Files, compressed_file_set, compression_processes)
    compress_sam_list(input_file_reader.controlRead1Files, compressed_file_set, compression_processes)
    compress_sam_list(input_file_reader.controlRead2Files, compressed_file_set, compression_processes)
    if len(compression_processes) > 0:
        print >> sys.stderr, "start to filter and compress sam files..."
        while utils.count_complete_process(compression_processes) < len(compression_processes):
            time.sleep(5)
        for proc, filename in compression_processes:
            status = proc.poll()
            if status == 0:
                print >> sys.stderr, '"%s" is filtered with only primary alignments and compressed to bam file by ' \
                                     'samtools.' % str(filename)
                subprocess.call(["rm", filename])
                print >> sys.stderr, '"%s" is removed.' % str(filename)
            else:
                print >> sys.stderr, 'error when compressing file "%s"' % str(filename)
        compression_end_time = time.time()
        print >> sys.stderr, 'cost %s minutes to compress all sam files.' % str(
            round((compression_end_time - compression_start_time) / 60.0, 2))


def compress_sam_list(filename_list, compressed_set, processes):
    """
    Compress a list of file and keep track the compression processes
    :param filename_list: list of file names
    :param compressed_set: set of files already compresses
    :param processes: list of compression processes
    :return:
    """
    for i in range(len(filename_list)):
        filename = filename_list[i]
        name = os.path.splitext(filename)[0]
        extension = os.path.splitext(filename)[-1].lower()
        if extension == ".sam":
            bam_name = name + ".bam"
            if filename not in compressed_set:
                processes.append(run_samtools(filename, bam_name))
                compressed_set.add(filename)
            filename_list[i] = bam_name


def run_samtools(sam_filename, bam_filename):
    """
    Run samtools to filter and compress into bam files
    :param sam_filename: name of the input sam file
    :param bam_filename: name of the output bam file
    :return: a samtools process to keep track
    """
    DEVNULL = open(os.devnull, 'w')
    f = open(bam_filename, 'w')
    p1_args = ["samtools", "view", "-S", "-F", "2048", "-h", sam_filename]
    p1 = subprocess.Popen(p1_args, stdout=subprocess.PIPE, stderr=DEVNULL)
    p2_args = ["samtools", "view", "-bS", "-"]
    p2 = subprocess.Popen(p2_args, stdin=p1.stdout, stdout=f, stderr=DEVNULL)
    p1.stdout.close()
    return p2, sam_filename
