#!/usr/bin/env python

"""
Generate matrices or vectors from Hi-C or control alignments
"""

import alabtools
import os
import datetime
import sys


class HicMatrix():
    """
    Create Hi-C contact matrix
    """
    def __init__(self, genome, resolution):
        self.genome = genome
        self.resolution = resolution
        info_file_path = os.path.join(os.path.dirname(alabtools.__file__), "genomes/" + genome + ".info")
        info_reader = open(info_file_path, 'r')
        info_content = info_reader.read().strip()
        self.chr_array = [x.split()[0] for x in info_content.split("\n")]
        self.length_array = [int(x.split()[1]) for x in info_content.split("\n")]
        start_idx_array = [0]
        for i in range(len(self.chr_array) - 1):
            start_idx_array.append(start_idx_array[i] + (self.length_array[i] - 1)/resolution + 1)
        self.chr_start_idx_dict = dict(zip(self.chr_array, start_idx_array))
        # construct the zero triangle contact matrix
        matrix_dimension = self.chr_start_idx_dict[self.chr_array[-1]] + (self.length_array[-1] - 1)/resolution + 1
        self.matrix = [[0] * (matrix_dimension - row) for row in range(matrix_dimension)]

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
        output_adj.write("#CHROM1\tSTART1\tEND1\tCHROM2\tSTART2\tEND2\tVALUE\n")
        # write file body
        for m in range(len(self.chr_array)):
            chr1 = self.chr_array[m]
            start_idx1 = self.chr_start_idx_dict[chr1]
            for i in range(start_idx1, len(self.matrix)):
                for n in range(m, len(self.chr_array)):
                    chr2 = self.chr_array[n]
                    start_idx2 = self.chr_start_idx_dict[chr2]
                    end_idx2 = start_idx2 + (self.length_array[n] - 1)/self.resolution + 1
                    for j in range(max(i, start_idx2), end_idx2):
                        contact_value = self.matrix[i][j]
                        if contact_value == 0:
                            pass
                        else:
                            chr1_idx = i - start_idx1
                            chr1_start = chr1_idx * self.resolution
                            chr1_end = min((chr1_idx + 1) * self.resolution, self.length_array[m])
                            chr2_idx = j - start_idx2
                            chr2_start = chr2_idx * self.resolution
                            chr2_end = min((chr2_idx + 1) * self.resolution, self.length_array[n])
                            output_adj.write("%s\t%d\t%d\t%s\t%d\t%d\t%d\n" % (chr1, chr1_start, chr1_end, chr2, chr2_start, chr2_end, contact_value))
        output_adj.close()


class CtlVector():
    """
    Create control bias vector
    """
    def __init__(self, genome, resolution):
        self.genome = genome
        self.resolution = resolution
        info_file_path = os.path.join(os.path.dirname(alabtools.__file__), "genomes/" + genome + ".info")
        info_reader = open(info_file_path, 'r')
        info_content = info_reader.read().strip()
        self.chr_array = [x.split()[0] for x in info_content.split("\n")]
        self.length_array = [int(x.split()[1]) for x in info_content.split("\n")]
        start_idx_array = [0]
        for i in range(len(self.chr_array) - 1):
            start_idx_array.append(start_idx_array[i] + (self.length_array[i] - 1)/resolution + 1)
        self.chr_start_idx_dict = dict(zip(self.chr_array, start_idx_array))
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
        output_bed.write("#CHROM\tSTART\tEND\tVALUE\n")
        # write file body
        for m in range(len(self.chr_array)):
            chrom = self.chr_array[m]
            start_idx = self.chr_start_idx_dict[chrom]
            end_idx = start_idx + (self.length_array[m] - 1)/self.resolution + 1
            for i in range(start_idx, end_idx):
                chr_idx = i - start_idx
                chr_start = chr_idx * self.resolution
                chr_end = min((chr_idx + 1) * self.resolution, self.length_array[m])
                output_bed.write("%s\t%d\t%d\t%d\n" % (chrom, chr_start, chr_end, self.vector[i]))
        output_bed.close()
