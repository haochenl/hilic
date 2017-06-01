#!/usr/bin/env python

"""
Hi-C data normalization workflow by experimental control
"""

import argparse
import os
import sys


__author__ = 'Haochen Li'
__version__ = 'v0.1.1'


class ProgramArguments():
    """
    The program command line parser
    """
    def __init__(self, description, version):
        self.parser = argparse.ArgumentParser(description=description, version=version)

    def parse(self, args=None):
        self.parser.add_argument("-c", "--chunk",
                                 help="the fastq file chunk size (number of lines, default: 90,000,000) must be times of 4",
                                 action="store", default=90000000, dest="chunk", metavar="INT", type=int)
        required_arguments = self.parser.add_argument_group("required arguments")
        required_arguments.add_argument("-r1", "--read1", help="path to the read1 fastq file",
                                        dest="read1", metavar="<PATH_TO_FILE>", required=True)
        required_arguments.add_argument("-r2", "--read2", help="path to the read2 fastq file",
                                        dest="read2", metavar="<PATH_TO_FILE>", required=True)
        required_arguments.add_argument("-o", "--output", help="the output directory prefix of the split fastq files",
                                        action="store", dest="outputPrefix", metavar="<PATH_TO_PREFIX>", required=True)
        return self.parser.parse_args(args)


def create_dir_and_fastqs(prefix, splits, read1_input, read2_input):
    directory = "%s_%d" % (prefix, splits)
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    out1 = open(os.path.join(directory, read1_input.split('.')[0] + ".%d.fastq" % splits), 'w')
    out2 = open(os.path.join(directory, read2_input.split('.')[0] + ".%d.fastq" % splits), 'w')
    return out1, out2


if __name__ == '__main__':
    ## parse command line arguments
    args = ProgramArguments(__doc__, __version__).parse()
    ## make split fastq directory and write paired fastq files
    chunk_reads = 0
    chunk_size = args.chunk / 4
    splits = 1
    read1 = open(args.read1, 'r')
    read2 = open(args.read2, 'r')
    r1_line = read1.readline()
    r2_line = read2.readline()
    out1, out2 = create_dir_and_fastqs(args.outputPrefix, splits, args.read1, args.read2)
    while r1_line and r2_line:
        if chunk_reads < chunk_size:
            if r1_line.split()[0][1:] == r2_line.split()[0][1:]:
                for i in range(4):
                    out1.write(r1_line)
                    out2.write(r2_line)
                    r1_line = read1.readline()
                    r2_line = read2.readline()
                chunk_reads += 1
            else:
                print >> sys.stderr, 'unmatched sequence query names between the paired fastq files.'
                sys.exit(1)
        else:
            out1.close()
            out2.close()
            splits += 1
            chunk_reads = 0
            out1, out2 = create_dir_and_fastqs(args.outputPrefix, splits, args.read1, args.read2)
    read1.close()
    read2.close()
    out1.close()
    out2.close()




