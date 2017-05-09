#!/usr/bin/env python

"""
Generate matrices or vectors from Hi-C or control alignments
"""

import alabtools
import alabtools.plots
import os
import datetime
import sys
import numpy as np


class HicMatrix():
    """
    Create Hi-C contact matrix
    """
    def __init__(self, genome, resolution, chrom_filter):
        self.genome = genome
        self.resolution = resolution
        info_file_path = os.path.join(os.path.dirname(alabtools.__file__), "genomes/" + genome + ".info")
        info_reader = open(info_file_path, 'r')
        info_content = info_reader.read().strip()
        self.chr_array = []
        self.length_array = []
        for line in info_content.split("\n"):
            fields = line.split()
            if fields[0] in chrom_filter:
                pass
            else:
                self.chr_array.append(fields[0])
                self.length_array.append(int(fields[1]))
        start_idx_array = [0]
        for i in range(len(self.chr_array) - 1):
            start_idx_array.append(start_idx_array[i] + (self.length_array[i] - 1)/resolution + 1)
        self.chr_start_idx_dict = dict(zip(self.chr_array, start_idx_array))
        print >> sys.stderr, '[generate Hi-C contact matrix]'
        # construct the zero triangle contact matrix
        matrix_dimension = self.chr_start_idx_dict[self.chr_array[-1]] + (self.length_array[-1] - 1)/resolution + 1
        self.matrix = [[0] * (matrix_dimension - row) for row in range(matrix_dimension)]
        print >> sys.stderr, 'initialized Hi-C contact matrix.'

    def populate(self, read1, read2):
        if read1.reference_name in self.chr_start_idx_dict and read2.reference_name in self.chr_start_idx_dict:
            read1_idx = self.chr_start_idx_dict[read1.reference_name] + (read1.pos - 1)/self.resolution
            read2_idx = self.chr_start_idx_dict[read2.reference_name] + (read2.pos - 1)/self.resolution
            if read1_idx <= read2_idx:
                self.matrix[read1_idx][read2_idx - read1_idx] += 1
            else:
                self.matrix[read2_idx][read1_idx - read2_idx] += 1
        else:
            pass

    def write(self, output_filename):
        print >> sys.stderr, 'write Hi-C contact matrix into contact adjacency list with resolution of %d' % self.resolution
        output_adj = open(output_filename, 'w')
        # write file header
        output_adj.write("##fileformat=ADJ\n")
        now = datetime.datetime.now()
        output_adj.write("##filedate=%d-%d-%d\n" % (now.year, now.month, now.day))
        output_adj.write("##resolution=%d\n" % self.resolution)
        for i in range(len(self.chr_array)):
            output_adj.write("##contig=<ID=%s,length=%d,assembly=%s>\n" % (self.chr_array[i], self.length_array[i], self.genome))
        output_adj.write("#CHROM1\tSTART1\tEND1\tINDEX1\tCHROM2\tSTART2\tEND2\tINDEX2\tVALUE\n")
        # write file body
        for m in range(len(self.chr_array)):
            chr1 = self.chr_array[m]
            start_idx1 = self.chr_start_idx_dict[chr1]
            end_idx1 = start_idx1 + (self.length_array[m] - 1)/self.resolution + 1
            for i in range(start_idx1, end_idx1):
                for n in range(m, len(self.chr_array)):
                    chr2 = self.chr_array[n]
                    start_idx2 = self.chr_start_idx_dict[chr2]
                    end_idx2 = start_idx2 + (self.length_array[n] - 1)/self.resolution + 1
                    for j in range(max(i, start_idx2), end_idx2):
                        contact_value = self.matrix[i][j-i]
                        if contact_value == 0:
                            pass
                        else:
                            chr1_idx = i - start_idx1
                            chr1_start = chr1_idx * self.resolution
                            chr1_end = min((chr1_idx + 1) * self.resolution, self.length_array[m])
                            chr2_idx = j - start_idx2
                            chr2_start = chr2_idx * self.resolution
                            chr2_end = min((chr2_idx + 1) * self.resolution, self.length_array[n])
                            output_adj.write("%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n" % (chr1, chr1_start, chr1_end, i, chr2, chr2_start, chr2_end, j, contact_value))
        output_adj.close()


class CtlVector():
    """
    Create control bias vector
    """
    def __init__(self, genome, resolution, chrom_filter):
        self.genome = genome
        self.resolution = resolution
        info_file_path = os.path.join(os.path.dirname(alabtools.__file__), "genomes/" + genome + ".info")
        info_reader = open(info_file_path, 'r')
        info_content = info_reader.read().strip()
        self.chr_array = []
        self.length_array = []
        for line in info_content.split("\n"):
            fields = line.split()
            if fields[0] in chrom_filter:
                pass
            else:
                self.chr_array.append(fields[0])
                self.length_array.append(int(fields[1]))
        start_idx_array = [0]
        for i in range(len(self.chr_array) - 1):
            start_idx_array.append(start_idx_array[i] + (self.length_array[i] - 1)/resolution + 1)
        self.chr_start_idx_dict = dict(zip(self.chr_array, start_idx_array))
        print >> sys.stderr, '[generate control bias vector]'
        # construct the zero control vector
        vector_dimension = self.chr_start_idx_dict[self.chr_array[-1]] + (self.length_array[-1] - 1)/resolution + 1
        self.vector = [0] * vector_dimension

    def populate(self, read):
        if not read.is_unmapped and read.reference_name in self.chr_start_idx_dict:
            read_idx = self.chr_start_idx_dict[read.reference_name] + (read.pos - 1)/self.resolution
            self.vector[read_idx] += 1
        else:
            pass

    def write(self, output_filename):
        print >> sys.stderr, 'write control vector into bed file with resolution of %d' % self.resolution
        output_bed = open(output_filename, 'w')
        # write file header
        output_bed.write("##fileformat=BED\n")
        now = datetime.datetime.now()
        output_bed.write("##filedate=%d-%d-%d\n" % (now.year, now.month, now.day))
        output_bed.write("##resolution=%d\n" % self.resolution)
        for i in range(len(self.chr_array)):
            output_bed.write("##contig=<ID=%s,length=%d,assembly=%s>\n" % (self.chr_array[i], self.length_array[i], self.genome))
        output_bed.write("#CHROM\tSTART\tEND\tINDEX\tVALUE\n")
        # write file body
        for m in range(len(self.chr_array)):
            chrom = self.chr_array[m]
            start_idx = self.chr_start_idx_dict[chrom]
            end_idx = start_idx + (self.length_array[m] - 1)/self.resolution + 1
            for i in range(start_idx, end_idx):
                chr_idx = i - start_idx
                chr_start = chr_idx * self.resolution
                chr_end = min((chr_idx + 1) * self.resolution, self.length_array[m])
                output_bed.write("%s\t%d\t%d\t%d\t%d\n" % (chrom, chr_start, chr_end, i, self.vector[i]))
        output_bed.close()


class MatrixNorm():
    def __init__(self, adj_filename, bed_filename):
        self.adj_filename = adj_filename
        self.bed_filename = bed_filename
        self.bias_vector = None
        self.contact_matrix = None

    def build_sss_matrix(self, dim):
        obj = []
        row = []
        col = []
        with open(self.adj_filename, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    pass
                else:
                    fields = line.strip().split()
                    row.append(int(fields[3]))
                    col.append(int(fields[7]))
                    obj.append(int(fields[8]))
        matrix_input = (obj, (row, col))
        return alabtools.matrix.sss_matrix(matrix_input, shape=(dim, dim))

    def build_bias_vector(self, cutoff=10):
        bias = []
        with open(self.bed_filename, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    pass
                else:
                    fields = line.strip().split()
                    value = int(fields[4])
                    if value < cutoff:
                        bias.append(0)
                    else:
                        bias.append(1000.0/value)
        self.bias_vector = np.asarray(bias).astype(np.float32)

    def create_raw_matrix(self, genome, resolution):
        contact_matrix = alabtools.api.Contactmatrix(genome=genome, resolution=resolution)
        self.build_bias_vector()
        contact_matrix.matrix = self.build_sss_matrix(len(self.bias_vector))
        self.contact_matrix = contact_matrix

    def ctlnorm(self):
        self.contact_matrix.matrix.norm(self.bias_vector)

    def krnorm(self):
        self.contact_matrix.maskLowCoverage()
        self.contact_matrix.krnorm()


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
