#!/usr/bin/env python

"""
Generate matrices or vectors from Hi-C or control alignments
"""

import alabtools
import os


class HicMatrix():
    """
    Create Hi-C contact matrix
    """
    def __init__(self, genome, resolution):
        self.resolution = resolution
        info_file_path = os.path.join(os.path.dirname(alabtools.__file__), "genomes/" + genome + ".info")
        info_content = info_file_path.read()
        chr_array = [x.split()[0] for x in info_content.split("\n")]
        length_array = [int(x.split()[1]) for x in info_content.split("\n")]
        start_idx_array = [0]
        for i in range(len(chr_array) - 1):
            start_idx_array.append(start_idx_array[i] + (length_array[i] - 1)/resolution + 1)
        self.chr_start_idx_dict = dict(zip(chr_array, start_idx_array))
        # construct the zero triangle contact matrix
        matrix_dimension = self.chr_start_idx_dict[chr_array[-1]] + (length_array[-1] - 1)/resolution + 1
        self.matrix = [[0] * (matrix_dimension - row) for row in range(matrix_dimension)]

    def populate(self, read1, read2):
        read1_idx = self.chr_start_idx_dict[read1.reference_name] + (read1.pos - 1)/self.resolution + 1
        read2_idx = self.chr_start_idx_dict[read2.reference_name] + (read2.pos - 1)/self.resolution + 1
        if read1_idx <= read2_idx:
            self.matrix[read1_idx][read2_idx - read1_idx] += 1
        else:
            self.matrix[read2_idx][read1_idx - read2_idx] += 1

class CtlVector():
    """
    Create control bias vector
    """
    def __init__(self, genome, resolution):
        self.resolution = resolution
        info_file_path = os.path.join(os.path.dirname(alabtools.__file__), "genomes/" + genome + ".info")
        info_content = info_file_path.read()
        chr_array = [x.split()[0] for x in info_content.split("\n")]
        length_array = [int(x.split()[1]) for x in info_content.split("\n")]
        start_idx_array = [0]
        for i in range(len(chr_array) - 1):
            start_idx_array.append(start_idx_array[i] + (length_array[i] - 1)/resolution + 1)
        self.chr_start_idx_dict = dict(zip(chr_array, start_idx_array))
        # construct the zero control vector
        vector_dimension = self.chr_start_idx_dict[chr_array[-1]] + (length_array[-1] - 1)/resolution + 1
        self.vector = [0] * vector_dimension

    def populate(self, read):
        if not read.is_unmapped:
            read_idx = self.chr_start_idx_dict[read.reference_name] + (read1.pos - 1)/self.resolution + 1
            self.vector[read_idx] += 1
        else:
            pass

