#!/usr/bin/env python

"""
Hi-C data normalization workflow by experimental control
"""

import argparse
import matrix
import sys
import plots


class ProgramArguments():
    """
    The program command line parser
    """
    def __init__(self, description, version):
        self.parser = argparse.ArgumentParser(description=description, version=version)

    def parse(self, args=None):
        required_arguments = self.parser.add_argument_group("required arguments")
        required_arguments.add_argument("-a", "--adj", help="path to the .adj file",
                                        dest="adj", metavar="<PATH_TO_FILE>", required=True)
        required_arguments.add_argument("-b", "--bed", help="path to the .bed file",
                                        dest="bed", metavar="<PATH_TO_FILE>", required=True)
        required_arguments.add_argument("-g", "--genome", help="reference genome of the input data", action="store",
                                        dest="genome", choices=["hg19", "hg38", "mm9", "mm10"],
                                        metavar="<hg19 OR hg38 OR mm9 OR mm10>",
                                        required=True)
        required_arguments.add_argument("-r", "--resolution",
                                        help="base pair resolution of the output contact matrix/matrices",
                                        action="store", dest="resolution", metavar="INT", type=int, required=True)
        required_arguments.add_argument("-o", "--output", help="the output prefix of the bam files",
                                        action="store", dest="outputPrefix", metavar="<PATH TO PREFIX>", required=True)
        return self.parser.parse_args(args)


if __name__ == '__main__':
    ## parse command line arguments
    args = ProgramArguments(__doc__, __version__).parse()
    ## generate Hi-C contact matrix
    print >> sys.stderr, '[load raw Hi-C matrix for resolution: %d]' % args,resolution
    raw_matrix = matrix.MatrixNorm(args.adj, args.bed)
    raw_matrix.create_raw_matrix(args.genome, args.resolution)
    raw_matrix.contact_matrix.save(args.outputPrefix + "_observed_%d" % args.resolution)
    print >> sys.stderr, 'plot heatmaps for individual chromosomes.'
    plots.plot_heatmaps(raw_matrix.contact_matrix, "heatmap_observed_%d" % args.resolution, "Observed")
    ## do Hi-C res normalization
    print >> sys.stderr, '[matrix normalization with Hi-C control bias file for resolution: %d]' % args.resolution
    hic_ctlnorm = copy.deepcopy(raw_matrix)
    hic_ctlnorm.ctlnorm()
    hic_ctlnorm.contact_matrix.save(args.outputPrefix + "_res_norm_%d" % args.resolution)
    print >> sys.stderr, 'plot heatmaps for individual chromosomes.'
    plots.plot_heatmaps(hic_ctlnorm.contact_matrix, "heatmap_res_norm_%d" % args.resolution, "RESnorm")
    ## do kr normalization
    print >> sys.stderr, '[matrix KR normalization for resolution: %d]' % args.resolution
    hic_krnorm = copy.deepcopy(raw_matrix)
    hic_krnorm.krnorm()
    hic_krnorm.contact_matrix.save(args.outputPrefix + "_kr_norm_%d" % args.resolution)
    print >> sys.stderr, 'plot heatmaps for individual chromosomes.'
    plots.plot_heatmaps(hic_krnorm.contact_matrix, "heatmap_kr_norm_%d" % res, "KRnorm")
    print >> sys.stderr, '[plot diff heatmaps between RES norm and KR norm]'
    plots.plot_diff_heatmaps("Diff", hic_ctlnorm.contact_matrix, "RES", hic_krnorm.contact_matrix, "KR")

