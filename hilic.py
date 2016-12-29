#!/usr/bin/env python

"""
Hi-C data normalization workflow by experimental input.
"""

import argparse

__author__ = 'H.L.'
__version__ = 'v0.0.1'

class ProgramArguments():
    def __init__(self, description, version):
        self.parser = argparse.ArgumentParser(description=description, version=version)

    def parse(self, args=None):
        self.parser.add_argument("-q", "--mapq", help="the MAPQ score filter threshold (default: lower than 30)", action="store", default=30, dest="mapq", metavar="INT", type=int)
        self.parser.add_argument("-o", "--observed", help="output observed contact matrix/matrices simultanously", action="store_true", dest="observed")
        self.parser.add_argument("-i", "--ice", help="output ICE normalized contact matrix/matrices simultanously", action="store_true", dest="ice")
        self.parser.add_argument("-l", "--log", help="output the program log file", action="store_true", dest="log")
        required_arguments = self.parser.add_argument_group("required arguments")
        required_arguments.add_argument("-c", "--config", help="path to the input configuration file", action="store", dest="config", metavar="<PATH_TO_FILE>", required=True)
        required_arguments.add_argument("-r", "--reference", help="reference genome of the input data", action="store", dest="reference", choices=["hg19", "hg38"], metavar="<hg19 OR hg38>", required=True)
        required_arguments.add_argument("-b", "--bin", help="base pair resolution(s) of the output contact matrix/matrices", action="store", nargs="+", dest="bin", metavar="INT", type=int, required=True)
        return self.parser.parse_args(args)

if __name__ == '__main__':
    args = ProgramArguments(__doc__, __version__).parse()
    print args.log
    print args.mapq
    print args.observed
    print args.ice
    print args.config
    print args.reference
    print args.bin
