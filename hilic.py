#!/usr/bin/env python

"""
Hi-C data normalization workflow by experimental control ([Hi]-C [Li]'s [C]orrection: HiLiC)
"""

import argparse
import sys
import re
import os
import linuxUtils

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
        self.parser.add_argument("-l", "--log", help="output the program log file", action="store_true", dest="log")
        required_arguments = self.parser.add_argument_group("required arguments")
        required_arguments.add_argument("-c", "--config", help="path to the input configuration file", action="store",
                                        dest="config", metavar="<PATH_TO_FILE>", required=True)
        required_arguments.add_argument("-g", "--genome", help="reference genome of the input data", action="store",
                                        dest="genome", choices=["hg19", "hg38"], metavar="<hg19 OR hg38>",
                                        required=True)
        required_arguments.add_argument("-r", "--resolution",
                                        help="base pair resolution(s) of the output contact matrix/matrices",
                                        action="store", nargs="+", dest="resolution", metavar="INT", type=int,
                                        required=True)

        return self.parser.parse_args(args)


class ConfigFileReader():

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
        self.hicRead1FastqFiles = []
        self.hicRead2FastqFiles = []
        self.hicRead1BamFiles = []
        self.hicRead2BamFiles = []
        self.controlRead1FastqFiles = []
        self.controlRead2FastqFiles = []
        self.controlRead1BamFiles = []
        self.controlRead2BamFiles = []

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

        hic_read1_format = os.path.splitext(hic_read1_file_list[0])[-1].lower()
        hic_read2_format = os.path.splitext(hic_read2_file_list[0])[-1].lower()
        control_read1_format = os.path.splitext(control_read1_file_list[0])[-1].lower()
        control_read2_format = os.path.splitext(control_read2_file_list[0])[-1].lower()

        if hic_read1_format == ".fastq":
            self.hicRead1FastqFiles = hic_read1_file_list
        elif hic_read1_format == ".bam":
            self.hicRead1BamFiles = hic_read1_file_list
        else:
            print >> sys.stderr, 'unsupported read1 Hi-C input format in config file '

        if hic_read2_format == ".fastq":
            self.hicRead2FastqFiles = hic_read2_file_list
        elif hic_read2_format == ".bam":
            self.hicRead2BamFiles = hic_read2_file_list
        else:
            print >> sys.stderr, 'unsupported read2 Hi-C input format in config file '

        if control_read1_format == ".fastq":
            self.controlRead1FastqFiles = control_read1_file_list
        elif control_read1_format == ".bam":
            self.controlRead1BamFiles = control_read1_file_list
        else:
            print >> sys.stderr, 'unsupported read1 control input format in config file '

        if control_read2_format == ".fastq":
            self.controlRead2FastqFiles = control_read2_file_list
        elif control_read2_format == ".bam":
            self.controlRead2BamFiles = control_read2_file_list
        else:
            print >> sys.stderr, 'unsupported read2 control input format in config file '

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


def check_file_status(fileList):
    for file in fileList:
        if os.path.isdir(file):
            print >> sys.stderr, 'directory path not supported in the config file'
            sys.exit(1)
        else:
            if not os.path.exists(file):
                print >> sys.stderr, 'listed file does not exist: %s' % str(file)
                sys.exit(1)


if __name__ == '__main__':
    args = ProgramArguments(__doc__, __version__).parse()
    config_reader = ConfigFileReader(args.config)
    config_reader.parse()
