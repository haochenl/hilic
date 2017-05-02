#!/usr/bin/env python

"""
Handle Hi-C or control datasets to retrieve control reads
"""

import pysam
import os
import sys
import time
import regex


class PairReads():
    _siteDictionary = {"HindIII": "AGCTT"}

    def __init__(self, read1_filename, read2_filename, type='rb'):
        self.read1 = pysam.AlignmentFile(read1_filename, type)
        self.read2 = pysam.AlignmentFile(read2_filename, type)
        self.read1_filename = read1_filename
        self.read2_filename = read2_filename

    def hic_separate(self, cutoff, mapq, enzyme, output_prefix):
        self.read1.reset()
        self.read2.reset()
        itr1 = self.read1.fetch(until_eof=True)
        itr2 = self.read2.fetch(until_eof=True)
        read1_name = os.path.splitext(self.read1_filename)[0]
        read2_name = os.path.splitext(self.read2_filename)[0]
        hic1_output = pysam.AlignmentFile(output_prefix + ".hic1.bam", 'wb', template=self.read1)
        hic2_output = pysam.AlignmentFile(output_prefix + ".hic2.bam", 'wb', template=self.read2)
        # the control reads alignment output
        ctl_output = pysam.AlignmentFile(output_prefix + ".ctl.bam", 'wb', template=self.read1)
        # the religation reads alignment output
        rlg_output = pysam.AlignmentFile(output_prefix + ".rlg.bam", 'wb', template=self.read1)
        # the single end aligned reads alignment output
        sgl_output = pysam.AlignmentFile(output_prefix + ".sgl.bam", 'wb', template=self.read1)
        # the junk alignment output
        junk_output = pysam.AlignmentFile(output_prefix + ".jk.bam", 'wb', template=self.read1)
        total = 0
        hic_count = 0
        ctl_count = 0
        rlg_count = 0
        sgl_count = 0
        junk_count = 0
        print >> sys.stderr, '[process Hi-C alignment files to retrieve contact reads and control reads]'
        start_time = time.time()
        for r1 in itr1:
            r2 = itr2.next()
            total += 1
            if r1.query_name == r2.query_name:
                # write the single end reads output
                if is_unmapped_or_low_mapq(r1, mapq) != is_unmapped_or_low_mapq(r2, mapq):
                    sgl_count += 1
                    sgl_output.write(r1)
                    sgl_output.write(r2)
                # write the junk output
                if is_junk(r1, r2, mapq):
                    junk_count += 1
                    junk_output.write(r1)
                    junk_output.write(r2)
                # write the Hi-C output
                elif is_hic(r1, r2, cutoff, mapq):
                    hic_count += 1
                    hic1_output.write(r1)
                    hic2_output.write(r2)
                else:
                    enzyme_site = self._siteDictionary[enzyme]
                    # write the normalization control output
                    if is_ctl(r1, enzyme_site) or is_ctl(r2, enzyme_site):
                        ctl_count += 1
                        ctl_output.write(r1)
                        ctl_output.write(r2)
                    # write the re-ligation reads into output
                    elif is_unmapped_or_low_mapq(r1, mapq) == is_unmapped_or_low_mapq(r2, mapq):
                        rlg_count += 1
                        rlg_output.write(r1)
                        rlg_output.write(r2)
            else:
                print >> sys.stderr, 'unmatched headers between two reads -> truncated files'
                sys.exit(1)
        print >> sys.stderr, 'total number of read pairs processed: %s' % str(total)
        print >> sys.stderr, 'number of Hi-C contacts identified: %s' % str(hic_pair_count)
        print >> sys.stderr, 'number of single end read pairs identified: %s' % str(sgl_count)
        print >> sys.stderr, 'number of control read pairs identified: %s' % str(ctl_count)
        print >> sys.stderr, 'number of re-ligation read pairs identified: %s' % str(rlg_count)
        print >> sys.stderr, 'number of junk read pairs: %s' % str(junk_count)
        hic1_output.close()
        hic2_output.close()
        ctl_output.close()
        rlg_output.close()
        sgl_output.close()
        junk_output.close()
        end_time = time.time()
        print >> sys.stderr, 'cost %s minutes to process Hi-C bam files.' % str(
            round((end_time - start_time) / 60.0, 2))

    def control_separate(self, cutoff, mapq, enzyme, output_prefix):
        self.read1.reset()
        self.read2.reset()
        itr1 = self.read1.fetch(until_eof=True)
        itr2 = self.read2.fetch(until_eof=True)
        ctl_output = pysam.AlignmentFile(output_prefix + ".ctl.bam", 'wb', template=self.read1)
        misc_output = pysam.AlignmentFile(output_prefix + ".misc.bam", 'wb', template=self.read1)
        junk_output = pysam.AlignmentFile(output_prefix + ".jk.bam", 'wb', template=self.read1)
        total = 0
        ctl_count = 0
        misc_count = 0
        junk_count = 0
        print >> sys.stderr, '[process control alignment files to retrieve control reads]'
        start_time = time.time()
        for r1 in itr1:
            r2 = itr2.next()
            total += 1
            if r1.query_name == r2.query_name:
                # write into junk output if read pair meets Hi-C conditions
                if is_hic(r1, r2, cutoff, mapq):
                    junk_count += 1
                    junk_output.write(r1)
                    junk_output.write(r2)
                else:
                    enzyme_site = self._siteDictionary[enzyme]
                    # write the normalization control output
                    if is_ctl(r1, enzyme_site) or is_ctl(r2, enzyme_site):
                        ctl_count += 1
                        ctl_output.write(r1)
                        ctl_output.write(r2)
                    # write the miscellaneous control read to output
                    else:
                        misc_count += 1
                        misc_output.write(r1)
                        misc_output.write(r2)
            else:
                print >> sys.stderr, 'unmatched headers between two reads -> truncated files'
                sys.exit(1)
        print >> sys.stderr, 'total number of read pairs processed: %s' % str(total)
        print >> sys.stderr, 'number of control read pairs identified: %s' % str(ctl_count)
        print >> sys.stderr, 'number of miscellaneous read pairs identified: %s' % str(misc_count)
        print >> sys.stderr, 'number of junk read pairs: %s' % str(junk_count)
        ctl_output.close()
        misc_output.close()
        junk_output.close()
        end_time = time.time()
        print >> sys.stderr, 'cost %s minutes to process control bam files.' % str(
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


def is_ctl(read, site):
    return regex.match("(^%s){e<=1}" % site, read.query_sequence) is not None
