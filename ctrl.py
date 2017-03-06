#!/usr/bin/env python

"""
Handle Hi-C or control datasets to retrieve control reads
"""

import pysam
import os
import sys


class PairReads():
    def __init__(self, read1_filename, read2_filename, type='rb'):
        self.read1 = pysam.AlignmentFile(read1_filename, type)
        self.read2 = pysam.AlignmentFile(read2_filename, type)

    def hic_separate(self, cutoff, mapq):
        self.read1.reset()
        self.read2.reset()
        itr1 = self.read1.fetch(until_eof=True)
        itr2 = self.read2.fetch(until_eof=True)
        read1_name = os.path.splitext(read1_name)[0]
        read2_name = os.path.splitext(read2_name)[0]
        hic1_output = pysam.AlignmentFile(read1_name + ".hic1.bam", 'wb', template=self.read1)
        hic2_output = pysam.AlignmentFile(read2_name + ".hic2.bam", 'wb', template=self.read2)
        ctl1_output = pysam.AlignmentFile(read1_name + ".ctl1.bam", 'wb', template=self.read1)
        ctl2_output = pysam.AlignmentFile(read2_name + ".ctl2.bam", 'wb', template=self.read2)
        junk1_output = pysam.AlignmentFile(read1_name + ".jk1.bam", 'wb', template=self.read1)
        junk2_output = pysam.AlignmentFile(read2_name + ".jk2.bam", 'wb', template=self.read2)
        total = 0
        hic_count = 0
        ctl_count = 0
        junk_count = 0
        print >> sys.stderr, '[process Hi-C alignment files to retrieve contact reads and control reads]'
        for r1 in itr1:
            r2 = itr2.next()
            total += 1
            if r1.query_name == r2.query_name:
                if is_junk(r1, r2, mapq):
                    junk_count += 1
                    junk1_output.write(r1)
                    junk2_output.write(r2)
                elif is_hic(r1, r2, cutoff, mapq):
                    hic_count += 1
                    hic1_output.write(r1)
                    hic2_output.write(r2)
                else:
                    ctl_count += 1
                    ctl1_output.write(r1)
                    ctl2_output.write(r2)
            else:
                print >> sys.stderr, 'unmatched headers between two reads -> truncated files'
                sys.exit(1)
        print >> sys.stderr, 'total number of paired reads processed: %s' % str(total)
        print >> sys.stderr, 'number of Hi-C contacts identified: %s' % str(hic_count)
        print >> sys.stderr, 'number of control paired reads identified: %s' % str(ctl_count)
        print >> sys.stderr, 'number of junk paired reads: %s' % str(junk_count)
        hic1_output.close()
        hic2_output.close()
        ctl1_output.close()
        ctl2_output.close()
        junk1_output.close()
        junk2_output.close()


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
