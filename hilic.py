#!/usr/bin/env python

"""
Hi-C data normalization workflow by experimental control ([Hi]-C [Li]'s [C]orrection: HiLiC)
"""

import argparse
import sys
import re
import os
import linuxUtils
import time
import subprocess

__author__ = 'Haochen Li'
__version__ = 'v0.0.1'


class ProgramArguments():
    def __init__(self, description, version):
        self.parser = argparse.ArgumentParser(description=description, version=version)

    def parse(self, args=None):
        self.parser.add_argument("-q", "--mapq",
                                 help="the MAPQ score alignment filter threshold (default: lower than 30)",
                                 action="store", default=30, dest="mapq", metavar="INT", type=int)
        self.parser.add_argument("-o", "--observed", help="output observed contact matrix/matrices simultanously",
                                 action="store_true", dest="observed")
        self.parser.add_argument("-i", "--ice", help="output ICE normalized contact matrix/matrices simultanously",
                                 action="store_true", dest="ice")
        ##self.parser.add_argument("-l", "--log", help="output the program log file", action="store_true", dest="log")
        required_arguments = self.parser.add_argument_group("required arguments")
        required_arguments.add_argument("-c", "--config", help="path to the input configuration file", action="store",
                                        dest="config", metavar="<PATH_TO_FILE>", required=True)
        required_arguments.add_argument("-g", "--genome", help="reference genome of the input data", action="store",
                                        dest="genome", choices=["hg19", "hg38"], metavar="<hg19 OR hg38>",
                                        required=True)
        required_arguments.add_argument("-f", "--fasta", help="path to the bwa indexed reference genome fasta file", action="store",
                                        dest="fasta", metavar="<PATH_TO_FILE>", required=True)
        required_arguments.add_argument("-r", "--resolution",
                                        help="base pair resolution(s) of the output contact matrix/matrices",
                                        action="store", nargs="+", dest="resolution", metavar="INT", type=int,
                                        required=True)
        return self.parser.parse_args(args)


class InputFileReader():
    _hicRead1Key = "HicRead1"
    _hicRead2Key = "HicRead2"
    _controlRead1Key = "ControlRead1"
    _controlRead2Key = "ControlRead2"

    def __init__(self, configFilePath):
        try:
            self.__reader = open(configFilePath, 'r')
        except Exception, e:
            print >> sys.stderr, "Exception: %s" % str(e)
            sys.exit(1)
        self.configFilePath = configFilePath
        self.hicRead1Files = []
        self.hicRead2Files = []
        self.controlRead1Files = []
        self.controlRead2Files = []

    def parse(self):
        file_content = self.__reader.read()
        file_dictionary = self.build_file_dictionary(file_content)
        hic_read1_file_list = file_dictionary[self._hicRead1Key]
        hic_read2_file_list = file_dictionary[self._hicRead2Key]
        control_read1_file_list = file_dictionary[self._controlRead1Key]
        control_read2_file_list = file_dictionary[self._controlRead2Key]
        os.chdir(os.path.dirname(self.configFilePath))
        check_file_status(hic_read1_file_list)
        check_file_status(hic_read2_file_list)
        check_file_status(control_read1_file_list)
        check_file_status(control_read2_file_list)
        self.hicRead1Files = hic_read1_file_list
        self.hicRead2Files = hic_read2_file_list
        self.controlRead1Files = control_read1_file_list
        self.controlRead2Files = control_read2_file_list

    def build_file_dictionary(self, content):
        content_list = filter(None, [x.strip() for x in re.split('\[|\]', content)])
        if len(content_list) != 8:
            print >> sys.stderr, 'unexpected config file format!'
            sys.exit(1)
        file_dictionary = {self._hicRead1Key: [], self._hicRead2Key: [], self._controlRead1Key: [],
                           self._controlRead2Key: []}
        for key in file_dictionary:
            if key not in content_list:
                print >> sys.stderr, 'expected header "[%s]" missing in the config file' % str(key)
                sys.exit(1)
        for index in range(len(content_list)):
            if content_list[index] in file_dictionary:
                file_dictionary[content_list[index]] = content_list[index + 1].split("\n")
        return file_dictionary

    def process_input_files(self, genome_fasta):
        ## check if a fastq file is already aligned
        ## useful when HiC and control files are the same
        print >> sys.stderr, '[start processing input files]'
        aligned_file_set = set()
        alignment_processes = []
        alignment_start_time = time.time()
        align_input_list(self.hicRead1Files, genome_fasta,  aligned_file_set, alignment_processes)
        align_input_list(self.hicRead2Files, genome_fasta, aligned_file_set, alignment_processes)
        align_input_list(self.controlRead1Files, genome_fasta, aligned_file_set, alignment_processes)
        align_input_list(self.controlRead2Files, genome_fasta, aligned_file_set, alignment_processes)
        if len(alignment_processes) > 0:
            print >> sys.stderr, 'start to align fastq files...'
            while count_complete_process(alignment_processes) < len(alignment_processes):
                time.sleep(5)
            for proc, filename in alignment_processes:
                status = proc.poll()
                if status == 0:
                    print >> sys.stderr, 'file "%s" is aligned by "bwa mem -5".' % str(filename)
                else:
                    print >> sys.stderr, 'error when aligning file "%s"' % str(filename)
            alignment_end_time = time.time()
            print >> sys.stderr, 'cost %s minutes to align all fastq files.' % str(round((alignment_end_time - alignment_start_time)/60.0, 2))
        ## check if a sam file is already compressed
        ## useful when HiC and control files are the same
        compressed_file_set = set()
        compression_processes = []
        compression_start_time = time.time()
        compress_sam_list(self.hicRead1Files, compressed_file_set, compression_processes)
        compress_sam_list(self.hicRead2Files, compressed_file_set, compression_processes)
        compress_sam_list(self.controlRead1Files, compressed_file_set, compression_processes)
        compress_sam_list(self.controlRead2Files, compressed_file_set, compression_processes)
        if len(compression_processes) > 0:
            print >> sys.stderr, "start to filter and compress sam files..."
            while count_complete_process(compression_processes) < len(compression_processes):
                time.sleep(5)
            for proc, filename in compression_processes:
                status = proc.poll()
                if status == 0:
                    print >> sys.stderr, '"%s" is filtered with only primary alignments and compressed to bam file by samtools.' % str(filename)
                    subprocess.call(["rm", filename])
                    print >> sys.stderr, '"%s" is removed.' % str(filename)
                else:
                    print >> sys.stderr, 'error when compressing file "%s"' % str(filename)
            compression_end_time = time.time()
            print >> sys.stderr, 'cost %s minutes to compress all sam files.' % str(round((compression_end_time - compression_start_time)/60.0, 2))


def compress_sam_list(filename_list, compressed_set, processes):
    for i in range(len(filename_list)):
        filename = filename_list[i]
        name = os.path.splitext(filename)[0]
        extension = os.path.splitext(filename)[-1].lower()
        if extension == ".sam":
            bam_name = name + ".bam"
            if filename not in compressed_set:
                processes.append(linuxUtils.run_samtools(filename, bam_name))
                compressed_set.add(filename)
            filename_list[i] = bam_name

def align_input_list(filename_list, genome_fasta, aligned_set, processes):
    for i in range(len(filename_list)):
        filename = filename_list[i]
        name = os.path.splitext(filename)[0]
        extension = os.path.splitext(filename)[-1].lower()
        if extension == ".fastq":
            sam_name = name + ".sam"
            if filename not in aligned_set:
                processes.append(linuxUtils.run_bwa(filename, genome_fasta, sam_name))
                aligned_set.add(filename)
            filename_list[i] = sam_name
        elif extension == ".sam":
            print >> sys.stderr, 'only hilic processed bam alignment files are supported as input.'
            sys.exit(1)
        elif extension == ".bam":
            print >> sys.stderr, '"%s" input file must follow hilic processed bam file format.'
            pass
        else:
            print >> sys.stderr, 'input file format (%s) not supported' % str(extension)
            sys.exit(1)

def check_file_status(filename_list):
    for filename in filename_list:
        if os.path.isdir(filename):
            print >> sys.stderr, 'file name is a directory not a regular file: %s' % str(filename)
            print >> sys.stderr, 'directory path not supported in the config file'
            sys.exit(1)
        else:
            if not os.path.exists(filename):
                print >> sys.stderr, 'listed file does not exist: %s' % str(filename)
                sys.exit(1)

def count_complete_process(processes):
    count = 0
    for proc, filename in processes:
        status = proc.poll()
        if status == None:
            continue
        else:
            count += 1
    return count


if __name__ == '__main__':
    args = ProgramArguments(__doc__, __version__).parse()
    input_reader = InputFileReader(args.config)
    input_reader.parse()
    input_reader.process_input_files(args.fasta)
