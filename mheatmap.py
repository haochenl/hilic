#!/usr/bin/env python

"""
Hi-C data normalization workflow by experimental control
"""

import argparse
import pypairix
import numpy as np
from alabtools.plots import *
from plots import plot_diff_heatmap
import re
import sys
import copy
import math


__author__ = 'Haochen Li'
__version__ = 'v0.0.1'


class ProgramArguments():
    """
    The program command line parser
    """
    def __init__(self, description, version):
        self.parser = argparse.ArgumentParser(description=description, version=version)

    def parse(self, args=None):
        self.parser.add_argument("-m", "--max", help="set the color maximum contact number in heatmap",
                                 action="store", default=None, dest="max", metavar="INT", type=int)
        self.parser.add_argument("-t", "--threshold", help="RES reads threshold for a bin",
                                 action="store", default=10, dest="threshold", metavar="INT", type=float)
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
                                        action="store", dest="resolution", metavar="INT", required=True, type=int)
        required_arguments.add_argument("-o", "--output", help="the output prefix of the heatmaps",
                                        action="store", dest="outputPrefix", metavar="<PATH TO PREFIX>", required=True)
        return self.parser.parse_args(args)


def plot(filename_prefix, pairs_filename, ctl_filename, query, resolution, clip_max, cutoff):
    query1 = query.split('|')[0]
    query2 = query.split('|')[1]
    query1_list = filter(None, re.split("[:\-]+", query1))
    query2_list = filter(None, re.split("[:\-]+", query2))
    chrom1 = query1_list[0]
    start1 = int(query1_list[1]) + 1
    end1 = int(query1_list[2])
    chrom2 = query2_list[0]
    start2 = int(query2_list[1]) + 1
    end2 = int(query2_list[2])

    rows = (end1 - start1) / resolution + 1
    columns = (end2 - start2) / resolution + 1
    shape = (rows, columns)
    matrix = np.zeros(shape)

    pairs = pypairix.open(pairs_filename)
    ctl = pypairix.open(ctl_filename)
    it_pairs = pairs.query2D(chrom1, start1, end1, chrom2, start2, end2)
    it1_ctl = ctl.query(chrom1, start1, end1)
    it2_ctl = ctl.query(chrom2, start2, end2)

    print >> sys.stderr, '[build contact matrix of the query region]'
    for pair in it_pairs:
        index1 = (int(pair[2]) - start1) / resolution
        index2 = (int(pair[4]) - start2) / resolution
        matrix[index1, index2] += 1
    raw_nonzero = matrix[np.nonzero(matrix)]
    raw_99qtl_scale = np.percentile(raw_nonzero, 95)##/np.max(raw_nonzero)
    ##raw_colormap = make_colormap([(1,1,1),(1,0,0),raw_99qtl_scale,(1,0,0),(0,0,0)], 'wrb')
    plotmatrix(filename_prefix + "_raw.pdf", matrix, cmap=red, clip_max=raw_99qtl_scale)

    ## get the matrix median before normalization
    matrix_kr = copy.copy(matrix)
    raw_median = np.median(raw_nonzero)

    print >> sys.stderr, '[normalize the contact matrix by RES control]'
    row_bias = np.zeros(rows)
    column_bias = np.zeros(columns)
    for ctl in it1_ctl:
        index = (int(ctl[2]) - start1) / resolution
        row_bias[index] += 1
    for ctl in it2_ctl:
        index = (int(ctl[2]) - start2) / resolution
        column_bias[index] += 1
    for row in range(rows):
        for column in range(row, columns):
            bias = get_bias(row_bias[row], column_bias[column], cutoff)
            if bias == 0.0:
                matrix[row, column] = 0.0
            else:
                matrix[row, column] = matrix[row, column]/bias
    res_nonzero = matrix[np.nonzero(matrix)]
    scale = raw_median/np.median(res_nonzero)
    matrix *= scale
    res_nonzero = matrix[np.nonzero(matrix)]
    res_99qtl_scale = np.percentile(res_nonzero, 95)##/np.max(res_nonzero)
    ##res_colormap = make_colormap([(1,1,1),(1,0,0),res_99qtl_scale,(1,0,0),(0,0,0)], 'wrb')
    plotmatrix(filename_prefix + "_res_norm.pdf", matrix, cmap=red, clip_max=res_99qtl_scale)

    print >> sys.stderr, '[normalize the contact matrix by KR bias vector]'
    number_kb = resolution / 1000
    row_kr_file = "/auto/cmb-08/fa/nhua/GSE63525/GM12878_combined/%dkb_resolution_intrachromosomal/%s/MAPQGE30/%s_%dkb.KRnorm" % (number_kb, chrom1, chrom1, number_kb)
    column_kr_file = "/auto/cmb-08/fa/nhua/GSE63525/GM12878_combined/%dkb_resolution_intrachromosomal/%s/MAPQGE30/%s_%dkb.KRnorm" % (number_kb, chrom2, chrom2, number_kb)
    with open(row_kr_file, 'r') as f:
        row_kr_lines = f.readlines()
    row_kr_bias = [float(e.strip()) for e in row_kr_lines]
    with open(column_kr_file, 'r') as f:
        column_kr_lines = f.readlines()
    column_kr_bias = [float(e.strip()) for e in column_kr_lines]
    row_start_bin = start1 / resolution
    column_start_bin = start2 / resolution
    for row in range(rows):
        for column in range(row, columns):
            row_index = row_start_bin + row
            column_index = column_start_bin + column
            if is_nan_or_zero(row_kr_bias[row_index]) or is_nan_or_zero(column_kr_bias[column_index]):
                matrix_kr[row, column] = 0.0
            else:
                matrix_kr[row, column] = matrix_kr[row, column]/(row_kr_bias[row_index] * column_kr_bias[column_index])
    scale = raw_median/np.median(matrix_kr[np.nonzero(matrix_kr)])
    matrix_kr *= scale
    kr_nonzero = matrix_kr[np.nonzero(matrix_kr)]
    kr_99qtl_scale = np.percentile(kr_nonzero, 95)##/np.max(kr_nonzero)
    ##kr_colormap = make_colormap([(1,1,1),(1,0,0),kr_99qtl_scale,(1,0,0),(0,0,0)], 'wrb')
    plotmatrix(filename_prefix + "_kr_norm.pdf", matrix_kr, cmap=red, clip_max=kr_99qtl_scale)

    print >> sys.stderr, '[plot diff heatmap between RES norm and KR norm]'
    plot_diff_heatmap(filename_prefix + "_diff.pdf", matrix, matrix_kr)

def get_bias(row_bias, column_bias, cutoff):
    if row_bias < cutoff or column_bias < cutoff:
        return 0.0
    else:
        return row_bias * column_bias
def is_nan_or_zero(number):
    return math.isnan(number) or number == 0.0


if __name__ == '__main__':
    ## parse command line arguments
    args = ProgramArguments(__doc__, __version__).parse()
    ## plot the observed and RES normalized matrices
    plot(args.outputPrefix, args.pairs, args.control, args.query, args.resolution, args.max, args.threshold)
