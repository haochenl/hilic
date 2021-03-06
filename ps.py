#!/usr/bin/env python

"""
Hi-C data normalization workflow by experimental control
"""

import argparse
import sys
import re
import os
import utils
import alignment
import samtools
import time
import ctrl
import matrix
import copy

__author__ = 'Haochen Li'
__version__ = 'v0.1.1'


class ProgramArguments():
    """
    The program command line parser
    """
    def __init__(self, description, version):
        self.parser = argparse.ArgumentParser(description=description, version=version)

    def parse(self, args=None):
        self.parser.add_argument("-q", "--mapq",
                                 help="the MAPQ score alignment filter threshold (default: lower than 30)",
                                 action="store", default=30, dest="mapq", metavar="INT", type=int)
        self.parser.add_argument("-l", "--length",
                                 help="the maximum insert length cutoff in bp (default: 1000)",
                                 action="store", default=1000, dest="len", metavar="INT", type=int)
        self.parser.add_argument("-p", "--processes",
                                 help="the number of additional threads to use (default: 4)",
                                 action="store", default=4, dest="processes", metavar="INT", type=int)
        required_arguments = self.parser.add_argument_group("required arguments")
        required_arguments.add_argument("-i", "--input", help="path to the paired bam file",
                                        dest="input", metavar="<PATH_TO_FILE>", required=True)
        required_arguments.add_argument("-t", "--type", help="type of the paired bam file",
                                        dest="type", choices=["hic", "ctl"], metavar="<hic OR ctl>", required=True)
        required_arguments.add_argument("-e", "--enzyme", help="name of the restriction enzyme (case sensitive e.g. HindIII)",
                                        action="store", dest="enzyme", metavar="STRING", required=True)
        required_arguments.add_argument("-o", "--output", help="the output prefix of the bam files",
                                        action="store", dest="outputPrefix", metavar="<PATH TO PREFIX>", required=True)
        return self.parser.parse_args(args)


class InputFileReader():
    """
    The program config file parser and processor
    """
    _enzymeKey = "Enzyme"
    _hicRead1Key = "HicRead1"
    _hicRead2Key = "HicRead2"
    _controlRead1Key = "ControlRead1"
    _controlRead2Key = "ControlRead2"
    _outputPrefixKey = "OutputPrefix"

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
        self.enzyme = None
        self.outputPrefix = None
        self.isControlSeparate = None

    def parse(self):
        """
        Parse the configurations and populate data file names
        """
        file_content = self.__reader.read()
        file_dictionary = self.build_file_dictionary(file_content)
        hic_read1_file_list = file_dictionary[self._hicRead1Key]
        hic_read2_file_list = file_dictionary[self._hicRead2Key]
        control_read1_file_list = file_dictionary[self._controlRead1Key]
        control_read2_file_list = file_dictionary[self._controlRead2Key]
        os.chdir(os.path.dirname(self.configFilePath))
        utils.check_file_status(hic_read1_file_list)
        utils.check_file_status(hic_read2_file_list)
        utils.check_file_status(control_read1_file_list)
        utils.check_file_status(control_read2_file_list)
        self.hicRead1Files = sorted(hic_read1_file_list)
        self.hicRead2Files = sorted(hic_read2_file_list)
        self.controlRead1Files = sorted(control_read1_file_list)
        self.controlRead2Files = sorted(control_read2_file_list)
        if self.hicRead1Files == self.controlRead1Files and self.hicRead2Files == self.controlRead2Files:
            self.isControlSeparate = False
        else:
            self.isControlSeparate = True
        if len(file_dictionary[self._enzymeKey]) != 1:
            print >> sys.stderr, 'hilic does not support multi-enzyme Hi-C experiment.'
            sys.exit(1)
        self.enzyme = file_dictionary[self._enzymeKey][0]
        if len(file_dictionary[self._outputPrefixKey]) != 1:
            print >> sys.stderr, 'only one output file prefix is allowed in the config file.'
            sys.exit(1)
        self.outputPrefix = file_dictionary[self._outputPrefixKey][0]

    def build_file_dictionary(self, content):
        """
        Parse the config file content build file name dictionary for different input file sections
        :param content: the config file content
        :return: file name dictionary for different sections
        """
        content_list = filter(None, [x.strip() for x in re.split('\[|\]', content)])
        if len(content_list) != 12:
            print >> sys.stderr, 'unexpected config file format!'
            sys.exit(1)
        file_dictionary = {self._hicRead1Key: [], self._hicRead2Key: [], self._controlRead1Key: [],
                           self._controlRead2Key: [], self._enzymeKey: [], self._outputPrefixKey: []}
        for key in file_dictionary:
            if key not in content_list:
                print >> sys.stderr, 'expected header "[%s]" missing in the config file' % str(key)
                sys.exit(1)
        for index in range(len(content_list)):
            if content_list[index] in file_dictionary:
                file_dictionary[content_list[index]] = content_list[index + 1].split("\n")
        return file_dictionary

    def process_input_files(self, genome_fasta):
        """
        Align the input files and compress the alignment output
        :param genome_fasta: the bwa aligner fasta file path
        """
        print >> sys.stderr, '[start processing input files]'
        alignment.align_from_reader(self, genome_fasta)

    def concat_input_files(self, genome):
        print >> sys.stderr, '[concatenate input files]'
        concat_start_time = time.time()
        hic_r1_name = "%s_hic_r1_%s.bam" % (self.outputPrefix, str(genome))
        hic_r2_name = "%s_hic_r2_%s.bam" % (self.outputPrefix, str(genome))
        if len(self.hicRead1Files) > 1:
            samtools.run_samtools_concatenation(self.hicRead1Files, hic_r1_name)
            self.hicRead1Files = [hic_r1_name]
        if len(self.hicRead2Files) > 1:
            samtools.run_samtools_concatenation(self.hicRead2Files, hic_r2_name)
            self.hicRead2Files = [hic_r2_name]
        if self.isControlSeparate:
            ctl_r1_name = "%s_ctl_r1_%s.bam" % (self.outputPrefix, str(genome))
            ctl_r2_name = "%s_ctl_r2_%s.bam" % (self.outputPrefix, str(genome))
            if len(self.controlRead1Files) > 1:
                samtools.run_samtools_concatenation(self.controlRead1Files, ctl_r1_name)
                self.controlRead1Files = [ctl_r1_name]
            if len(self.controlRead2Files) > 1:
                samtools.run_samtools_concatenation(self.controlRead2Files, ctl_r2_name)
                self.controlRead2Files = [ctl_r2_name]
        else:
            self.controlRead1Files = [hic_r1_name]
            self.controlRead2Files = [hic_r2_name]
        concat_end_time = time.time()
        print >> sys.stderr, 'cost %s minutes to combine all bam files.' % str(
            round((concat_end_time - concat_start_time) / 60.0, 2))


if __name__ == '__main__':
    ## parse command line arguments
    args = ProgramArguments(__doc__, __version__).parse()
    if args.input != "-":
        utils.check_file_status(args.input)
    ## if not build directly from separate bam files, process the configuration input files
    if args.type == "hic":
        ## separate Hi-C pairs into contacts and control
        hic_pair = ctrl.PairReads(args.input)
        hic_pair.hic_separate(args.len, args.mapq, args.enzyme, args.outputPrefix, args.processes)
    if args.type == "ctl":
        ## process control experiment files if exists
        ctl_pair = ctrl.PairReads(args.input)
        ctl_pair.control_separate(args.len, args.mapq, args.enzyme, args.outputPrefix, args.processes)
