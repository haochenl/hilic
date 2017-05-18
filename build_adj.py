#!/usr/bin/env python

"""
Hi-C data normalization workflow by experimental control
"""

import argparse
import matrix
import sys

__author__ = 'Haochen Li'
__version__ = 'v0.0.1'

class ProgramArguments():
    """
    The program command line parser
    """
    def __init__(self, description, version):
        self.parser = argparse.ArgumentParser(description=description, version=version)

    def parse(self, args=None):
        required_arguments = self.parser.add_argument_group("required arguments")
        required_arguments.add_argument("-p", "--pairs", help="path to the pairs file",
                                        dest="pairs", metavar="<PATH_TO_FILE>", required=True)
        required_arguments.add_argument("-c", "--control", help="path to the control file",
                                        dest="control", metavar="<PATH_TO_FILE>", required=True)
        required_arguments.add_argument("-r", "--resolution", help="output matrix resolution",
                                        action="store", dest="resolution", metavar="INT", required=True, type=int)
        required_arguments.add_argument("-g", "--genome", help="reference genome of the input data",
                                        action="store", choices=["hg19", "hg38", "mm9", "mm10"], metavar="<hg19 OR hg38 OR mm9 OR mm10>", required=True)
        required_arguments.add_argument("-o", "--output", help="the output prefix of the bam files",
                                        action="store", dest="outputPrefix", metavar="<PATH TO PREFIX>", required=True)
        return self.parser.parse_args(args)


if __name__ == '__main__':
    ## parse command line arguments
    args = ProgramArguments(__doc__, __version__).parse()
    ## chromosomes filtered from the genome info resource file
    filter_set = {"chrY", "chrM"}
    print >> sys.stderr, "[build observed contact adjacency list at %d resolution]" % args.resolution
    hic_matrix = matrix.HicMatrix(args.genome, args.resolution, filter_set)
    with open(args.pairs) as f:
        for line in f:
            hic_matrix.populate(line)
    hic_matrix.write(args.outputPrefix + "_hic_%d.adj" % args.resolution)
    print >> sys.stderr, "[build control bias vector bed file at %d resolution]" % args.resolution
    hic_bias_vector = matrix.CtlVector(args.genome, args.resolution, filter_set)
    with open(args.control) as f:
        for line in f:
            hic_bias_vector.populate(line)
    hic_bias_vector.write(args.outputPrefix + "_ctl_%d.bed" % args.resolution)


