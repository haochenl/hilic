#!/usr/bin/env python

"""
Hi-C data normalization workflow by experimental control
"""

import argparse
import pypairix
import numpy as np
import alabtools.plots as plt
import re
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
        required_arguments.add_argument("-p", "--pairs", help="path to the pairix indexed pair file",
                                        dest="pairs", metavar="<PATH_TO_FILE>", required=True)
        required_arguments.add_argument("-c", "--control", help="path to the pairix indexed control file",
                                        dest="control", metavar="<PATH_TO_FILE>", required=True)
        required_arguments.add_argument("-q", "--query", help="the pairix style query string",
                                        dest="query",  metavar="STRING", required=True)
        required_arguments.add_argument("-g", "--genome", help="the reference genome of interest", dest="genome",
                                        choices=["hg19", "hg38", "mm9", "mm10"], metavar="<hg19 OR hg38 OR mm9 OR mm10>", required=True)
        required_arguments.add_argument("-r", "--resolution", help="the resolution of output heatmaps",
                                        action="store", dest="resolution", metavar="INT", required=True)
        required_arguments.add_argument("-o", "--output", help="the output prefix of the heatmaps",
                                        action="store", dest="outputPrefix", metavar="<PATH TO PREFIX>", required=True)
        return self.parser.parse_args(args)


def plot(filename_prefix, pairs_filename, ctl_filename, query, resolution):
    query1 = query.split('|')[0]
    query2 = query.split('|')[1]
    query1_list = filter(None, re.split("[:\-]+", query1))
    query2_list = filter(None, re.split("[:\-]+", query2))
    chrom1 = query1_list[0]
    start1 = int(query1_list[1])
    end1 = int(query1_list[2])
    chrom2 = query2_list[0]
    start2 = int(query2_list[1])
    end2 = int(query2_list[2])

    rows = (end1 - start1) / resolution
    columns = (end2 - start2) / resolution
    shape = (rows, columns)
    matrix = np.zeros(shape)

    pairs = pypairix.open(pairs_filename)
    ctl = pypairix.open(ctl_filename)
    it_pairs = pairs.querys2D(query)
    it1_ctl = ctl.querys2D(query1)
    it2_ctl = ctl.querys2D(query2)

    print >> sys.stderr, '[build contact matrix of the query region]'
    for pair in it_pairs:
        index1 = (int(pair[2]) - start1 - 1) / resolution
        index2 = (int(pair[4]) - start2 - 1) / resolution
        matrix[index1, index2] += 1
    plt.plotmatrix(filename_prefix + "_raw.pdf", matrix)

    print >> sys.stderr, '[normalize the contact matrix by RES control]'
    row_bias = np.zeros(rows)
    column_bias = np.zeros(columns)
    for ctl in it1_ctl:
        index = (int(ctl[2] - start1 - 1)) / resolution
        row_bias[index] += 1
    for ctl in it2_ctl:
        index = (int(ctl[2] - start2 - 1)) / resolution
        column_bias[index] += 1
    for row in range(rows):
        for column in range(row, columns):
            matrix[row, column] = matrix[row, column]/(row_bias[row]*column_bias[column])
    plt.plotmatrix(filename_prefix + "_res_norm.pdf", matrix)


if __name__ == '__main__':
    ## parse command line arguments
    args = ProgramArguments(__doc__, __version__).parse()
    ## plot the observed and RES normalized matrices
    plot(args.outputPrefix, args.pairs, args.control, args.query, args.resolution)



