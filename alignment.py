#!/usr/bin/env python

"""
The hilic alignment utilities
"""

import os
import sys
import time
import subprocess
import utils


def align_from_reader(input_file_reader, genome_fasta):
    """
    Align the fastq files from hilic config file reader
    :param input_file_reader: the hilic config file reader
    :param genome_fasta: the bwa aligner fasta file path
    """
    ## check if a fastq file is already aligned
    ## useful when HiC and control files are the same
    aligned_file_set = set()
    alignment_processes = []
    alignment_start_time = time.time()
    align_input_list(input_file_reader.hicRead1Files, genome_fasta, aligned_file_set, alignment_processes)
    align_input_list(input_file_reader.hicRead2Files, genome_fasta, aligned_file_set, alignment_processes)
    align_input_list(input_file_reader.controlRead1Files, genome_fasta, aligned_file_set, alignment_processes)
    align_input_list(input_file_reader.controlRead2Files, genome_fasta, aligned_file_set, alignment_processes)
    if len(alignment_processes) > 0:
        print >> sys.stderr, 'start to align fastq files...'
        while utils.count_complete_process(alignment_processes) < len(alignment_processes):
            time.sleep(5)
        for proc, filename in alignment_processes:
            status = proc.poll()
            if status == 0:
                print >> sys.stderr, 'file "%s" is aligned by "bwa mem -5".' % str(filename)
            else:
                print >> sys.stderr, 'error when aligning file "%s"' % str(filename)
        alignment_end_time = time.time()
        print >> sys.stderr, 'cost %s minutes to align all fastq files.' % str(
            round((alignment_end_time - alignment_start_time) / 60.0, 2))


def align_input_list(filename_list, genome_fasta, aligned_set, processes):
    """
    Run bwa aligner given a list of file names and add processes to a list
    :param filename_list: list of file names
    :param genome_fasta: the bwa aligner fasta file path
    :param aligned_set: set of files that already aligned
    :param processes: list of processes to track
    """
    for i in range(len(filename_list)):
        filename = filename_list[i]
        name = os.path.splitext(filename)[0]
        extension = os.path.splitext(filename)[-1].lower()
        if extension == ".fastq":
            sam_name = name + ".sam"
            if filename not in aligned_set:
                processes.append(run_bwa(filename, genome_fasta, sam_name))
                aligned_set.add(filename)
            filename_list[i] = sam_name
        elif extension == ".sam":
            print >> sys.stderr, 'only hilic processed bam alignment files are supported as input.'
            sys.exit(1)
        elif extension == ".bam":
            print >> sys.stderr, '"%s" input file must follow hilic processed bam file format. can\'t guarantee ' \
                                 'correct results if not.' % str(filename)
            pass
        else:
            print >> sys.stderr, 'input file format (%s) not supported' % str(extension)
            sys.exit(1)


def run_bwa(fastq_filename, genome_fasta, sam_filename):
    """
    Run bwa aligner and return an alignment process to track
    :param fastq_filename: the input fastq file name
    :param genome_fasta: the bwa aligner fasta file path
    :param sam_filename: the output sam file name
    :return: a bwa alignment process
    """
    DEVNULL = open(os.devnull, 'w')
    f = open(sam_filename, 'w')
    args = ["bwa", "mem", "-5", genome_fasta, fastq_filename]
    process = subprocess.Popen(args, stdout=f, stderr=DEVNULL)
    return process, fastq_filename
