#!/usr/bin/env python

"""
Handle Hi-C or control datasets to retrieve control reads
"""

import pysam
import os
import sys
import time


class PairReads():
    def __init__(self, read1_filename, read2_filename, type='rb'):
        self.read1 = pysam.AlignmentFile(read1_filename, type)
        self.read2 = pysam.AlignmentFile(read2_filename, type)
        self.read1_filename = read1_filename
        self.read2_filename = read2_filename

    def hic_separate(self, cutoff, mapq):
        self.read1.reset()
        self.read2.reset()
        itr1 = self.read1.fetch(until_eof=True)
        itr2 = self.read2.fetch(until_eof=True)
        read1_name = os.path.splitext(self.read1_filename)[0]
        read2_name = os.path.splitext(self.read2_filename)[0]
        output_prefix = "hic_combined"
        hic1_output = pysam.AlignmentFile(read1_name + ".hic1.bam", 'wb', template=self.read1)
        hic2_output = pysam.AlignmentFile(read2_name + ".hic2.bam", 'wb', template=self.read2)
        ctl_output = pysam.AlignmentFile(output_prefix + ".ctl.bam", 'wb', template=self.read1)
        junk_output = pysam.AlignmentFile(output_prefix + ".jk.bam", 'wb', template=self.read1)
        total = 0
        hic_pair_count = 0
        ctl_count = 0
        junk_count = 0
        print >> sys.stderr, '[process Hi-C alignment files to retrieve contact reads and control reads]'
        start_time = time.time()
        for r1 in itr1:
            r2 = itr2.next()
            total += 2
            if r1.query_name == r2.query_name:
                if is_junk(r1, r2, mapq):
                    junk_count += 2
                    junk_output.write(r1)
                    junk_output.write(r2)
                elif is_hic(r1, r2, cutoff, mapq):
                    hic_pair_count += 1
                    hic1_output.write(r1)
                    hic2_output.write(r2)
                else:
                    if is_unmapped_or_low_mapq(r1, mapq):
                        junk_count += 1
                        junk_output.write(r1)
                    else:
                        ctl_count += 1
                        ctl_output.write(r1)
                    if is_unmapped_or_low_mapq(r2, mapq):
                        junk_count += 1
                        junk_output.write(r2)
                    else:
                        ctl_count += 1
                        ctl_output.write(r2)
            else:
                print >> sys.stderr, 'unmatched headers between two reads -> truncated files'
                sys.exit(1)
        print >> sys.stderr, 'total number of reads processed: %s' % str(total)
        print >> sys.stderr, 'number of Hi-C contacts identified: %s' % str(hic_pair_count)
        print >> sys.stderr, 'number of control reads identified: %s' % str(ctl_count)
        print >> sys.stderr, 'number of junk reads: %s' % str(junk_count)
        hic1_output.close()
        hic2_output.close()
        ctl_output.close()
        junk_output.close()
        end_time = time.time()
        print >> sys.stderr, 'cost %s minutes to process all Hi-C bam files.' % str(
                            round((end_time - start_time) / 60.0, 2))


def is_unmapped_or_low_mapq(read, mapq):
    return read.is_unmapped or read.mapping_quality < mapq


def is_junk(read1, read2, mapq):
    return is_unmapped_or_low_mapq(read1, mapq) and is_unmapped_or_low_mapq(read2, mapq)


def is_hic(read1, read2, cutoff, mapq):
    if is_unmapped_or_low_mapq(read1, mapq) or is_unmapped_or_low_mapq(read2, mapq):
        return False
    if read1.reference_id != read2.reference_id:
        return True
    else:
        if abs(read1.pos - read2.pos) > cutoff:
            return True
        else:
            return False
