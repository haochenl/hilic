#!/usr/bin/env python

"""
Hi-C Data normalization workflow by experimental input.
"""

import argparse

__author__ = 'H.L.'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-l", "--logger", help="print the program log file", action="store_true", dest="logger")
    parser.add_argument("-q", "--mapq", help="the MAPQ score filter threshold (default: 30)", action="store", default=30, dest="mapq", metavar="INT", type=int)
    required_arguments = parser.add_argument_group("required arguments")
    required_arguments.add_argument("-c", "--config", help="path to the input configuration file", dest="config", metavar="<PATH_TO_FILE>", required=True)
    args = parser.parse_args()
    print args.logger
    print args.mapq
    print args.config