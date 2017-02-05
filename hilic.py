#!/usr/bin/env python

"""
Hi-C data normalization workflow by experimental control ([Hi]-C [Li]'s [C]orrection: HiLiC)
"""

import argparse
import sys
import re
import os
import utils
import alignment
import samtools

__author__ = 'Haochen Li'
__version__ = 'v0.0.1'


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
        required_arguments.add_argument("-f", "--fasta", help="path to the bwa indexed reference genome fasta file",
                                        action="store",
                                        dest="fasta", metavar="<PATH_TO_FILE>", required=True)
        required_arguments.add_argument("-r", "--resolution",
                                        help="base pair resolution(s) of the output contact matrix/matrices",
                                        action="store", nargs="+", dest="resolution", metavar="INT", type=int,
                                        required=True)
        return self.parser.parse_args(args)


class InputFileReader():
    """
    The program config file parser and processor
    """
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
        self.hicRead1Files = hic_read1_file_list
        self.hicRead2Files = hic_read2_file_list
        self.controlRead1Files = control_read1_file_list
        self.controlRead2Files = control_read2_file_list

    def build_file_dictionary(self, content):
        """
        Parse the config file content build file name dictionary for different input file sections
        :param content: the config file content
        :return: file name dictionary for different sections
        """
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
        """
        Align the input files and compress the alignment output
        :param genome_fasta: the bwa aligner fasta file path
        """
        print >> sys.stderr, '[start processing input files]'
        alignment.align_from_reader(self, genome_fasta)
        samtools.compress_from_reader(self)


if __name__ == '__main__':
    args = ProgramArguments(__doc__, __version__).parse()
    input_reader = InputFileReader(args.config)
    input_reader.parse()
    input_reader.process_input_files(args.fasta)
