"""
The hilic plotting utilities
"""

import alabtools
import numpy as np
import alabtools.plots
import math
import os

def custom_log(a, b):
    if a == 0.0 or b == 0.0:
        return 0.0
    else:
        return math.log(a/b)

def plot_diff_heatmap(filename, matrix1, matrix2, title=None):
    """
    compute and plot log ratio difference between matrices applied with different methods
    :param matrix1: a 2D numpy array
    :param method1: a string represents the feature of matrix1
    :param matrix2: a 2D numpy array
    :param method2: a string represents the feature of matrix2
    """
    median1 = np.median(matrix1[np.nonzero(matrix1)])
    median2 = np.median(matrix2[np.nonzero(matrix2)])
    ## normlize matrix1 to make it have same median as matrix2
    matrix1 *= (median2/median1)
    ## calculate the log ratio diff matrix
    matrix_log = np.vectorize(custom_log)
    matrix_diff = matrix_log(matrix1, matrix2)
    diff_min = np.min(matrix_diff)
    diff_max = np.max(matrix_diff)
    scale = float(diff_min/(diff_min - diff_max))
    if diff_min + diff_max < 0:
        redness = float(-diff_max/diff_min)
        color = alabtools.plots.make_colormap([(0.0,0.0,1.0),(1.0,1.0,1.0),scale,(1.0,1.0,1.0),(1.0,1.0-redness,1.0-redness)], 'bwr')
        alabtools.plots.plotmatrix(filename, matrix_diff, cmap=color, title=title)
    else:
        blueness = float(-diff_min/diff_max)
        color = alabtools.plots.make_colormap([(1.0-blueness,1.0-blueness,1.0),(1.0,1.0,1.0),scale,(1.0,1.0,1.0),(1.0,0.0,0.0)], 'bwr')
        alabtools.plots.plotmatrix(filename, matrix_diff, cmap=color, title=title)


def plot_heatmaps(contact_matrix, figure_prefix, method):
    directory = os.path.join(os.getcwd(), figure_prefix)
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    stella = alabtools.plots.make_colormap([(0,0,0),(0,0,1),0.0,(0,0,1),(1,1,0),0.05,(1,1,0),(1,1,1)],'stella')
    for chrom in contact_matrix.genome.chroms:
        matrix = contact_matrix[chrom].matrix.toarray()
        alabtools.plots.plotmatrix(os.path.join(directory, figure_prefix + "_stella_%s.pdf" % chrom), matrix, cmap=stella, title="%s(%s)" % (chrom, method))
        nonzeros = matrix[np.nonzero(matrix)]
        scale_95qtl = np.percentile(nonzeros, 95)
        alabtools.plots.plotmatrix(os.path.join(directory, figure_prefix + "_redness_%s.pdf" % chrom), matrix, clip_max=scale_95qtl, cmap=alabtools.plots.red, title="%s(%s)" % (chrom, method))


def plot_diff_heatmaps(figure_prefix, contact_matrix1, method1, contact_matrix2, method2):
    directory = os.path.join(os.getcwd(), figure_prefix)
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    for chrom in contact_matrix1.genome.chroms:
        matrix1 = contact_matrix1[chrom].matrix.toarray()
        matrix2 = contact_matrix2[chrom].matrix.toarray()
        plot_diff_heatmap(os.path.join(directory, figure_prefix + "_%s.pdf" % chrom), matrix1, matrix2, title="%s: log(%s/%s)" % (chrom, method1, method2))
