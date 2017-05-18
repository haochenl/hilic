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

def plot_diff_heatmaps(filename, matrix1, matrix2, title=None):
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
        alabtools.plots.plotmatrix(filename, mat_diff, cmap=color, title=title)


def plot_heatmaps(contact_matrix, figure_prefix, method, clip_max=None):
    directory = os.path.join(os.getcwd(), figure_prefix)
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    stella = alabtools.plots.make_colormap([(0,0,0),(0,0,1),0.0,(0,0,1),(1,1,0),0.05,(1,1,0),(1,1,1)],'stella')
    bloody = alabtools.plots.make_colormap([(1,1,1),(1,0,0),0.5,(1,0,0),(0,0,0)],'bloody')
    for chrom in contact_matrix.genome.chroms:
        contact_matrix[chrom].plot(os.path.join(directory, figure_prefix + "_stella_%s.pdf" % chrom), cmap=stella, title="%s(%s)" % (chrom, method))
        contact_matrix[chrom].plot(os.path.join(directory, figure_prefix + "_bloody_%s.pdf" % chrom), clip_max=clip_max, cmap=bloody, title="%s(%s)" % (chrom, method))
