#!/usr/bin/env python

"""
Handle Hi-C or control datasets to retrieve control reads
"""

import pysam
import os
import sys
import time
import regex
import matrix

class PairReads():
    _siteDictionary = {"HindIII": "AGCTT", "DpnII": "GATC", "NcoI": "CATGG", "MboI": "GATC"}
    _junctionDictionary = {"HindIII": "AAGCTAGCTT", "DpnII": "GATCGATC", "NcoI": "CCATGCATGG", "MboI": "GATCGATC"}

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
        hic_output_prefix = output_prefix + "_hic"
        hic1_output = pysam.AlignmentFile(hic_output_prefix + ".hic1.bam", 'wb', template=self.read1)
        hic2_output = pysam.AlignmentFile(hic_output_prefix + ".hic2.bam", 'wb', template=self.read2)
        # the control reads alignment output
        ctl_output = pysam.AlignmentFile(hic_output_prefix + ".ctl.bam", 'wb', template=self.read1)
        # the religation reads alignment output
        rlg_output = pysam.AlignmentFile(hic_output_prefix + ".rlg.bam", 'wb', template=self.read1)
        # the single end aligned reads alignment output
        sgl_output = pysam.AlignmentFile(hic_output_prefix + ".sgl.bam", 'wb', template=self.read1)
        # the junk alignment output
        junk_output = pysam.AlignmentFile(hic_output_prefix + ".jk.bam", 'wb', template=self.read1)
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
                is_r1_invalid = is_unmapped_or_low_mapq(r1, mapq)
                is_r2_invalid = is_unmapped_or_low_mapq(r2, mapq)
                # write the single end reads output
                if is_r1_invalid != is_r2_invalid:
                    sgl_count += 1
                    if not is_r1_invalid:
                        sgl_output.write(r1)
                    if not is_r2_invalid:
                        sgl_output.write(r2)
                if is_r1_invalid and is_r2_invalid:
                    junk_count += 1
                    junk_output.write(r1)
                    junk_output.write(r2)
                # write the Hi-C output
                elif is_hic(r1, r2, cutoff, is_r1_invalid, is_r2_invalid):
                    hic_count += 1
                    hic1_output.write(r1)
                    hic2_output.write(r2)
                else:
                    enzyme_site = self._siteDictionary[enzyme]
                    ligation_junction = self._junctionDictionary[enzyme]
                    # write the normalization control output
                    if is_ctl(r1, r2, enzyme_site, ligation_junction, is_r1_invalid, is_r2_invalid):
                        ctl_count += 1
                        if not is_r1_invalid:
                            ctl_output.write(r1)
                        if not is_r2_invalid:
                            ctl_output.write(r2)
                    # write the re-ligation reads into output
                    elif is_rlg(r1, r2, is_r1_invalid, is_r2_invalid):
                        rlg_count += 1
                        rlg_output.write(r1)
                        rlg_output.write(r2)
                    else:
                        junk_count += 1
                        junk_output.write(r1)
                        junk_output.write(r2)
            else:
                print >> sys.stderr, 'unmatched headers between two reads -> truncated files'
                sys.exit(1)
        print >> sys.stderr, 'total number of read pairs processed: %s' % str(total)
        print >> sys.stderr, 'number of Hi-C contacts identified: %s' % str(hic_count)
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
        ctl_output_prefix = output_prefix + "_ctl"
        ctl_output = pysam.AlignmentFile(ctl_output_prefix + ".ctl.bam", 'wb', template=self.read1)
        misc_output = pysam.AlignmentFile(ctl_output_prefix + ".misc.bam", 'wb', template=self.read1)
        junk_output = pysam.AlignmentFile(ctl_output_prefix + ".jk.bam", 'wb', template=self.read1)
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
                is_r1_invalid = is_unmapped_or_low_mapq(r1, mapq)
                is_r2_invalid = is_unmapped_or_low_mapq(r2, mapq)
                # write into junk output if read pair meets Hi-C conditions
                if is_hic(r1, r2, cutoff, is_r1_invalid, is_r2_invalid):
                    junk_count += 1
                    junk_output.write(r1)
                    junk_output.write(r2)
                else:
                    enzyme_site = self._siteDictionary[enzyme]
                    ligation_junction = self._junctionDictionary[enzyme]
                    # write the normalization control output
                    if is_ctl(r1, r2, enzyme_site, ligation_junction, is_r1_invalid, is_r2_invalid):
                        ctl_count += 1
                        if not is_r1_invalid:
                            ctl_output.write(r1)
                        if not is_r2_invalid:
                            ctl_output.write(r2)
                    # write the miscellaneous control read to output
                    else:
                        misc_count += 1
                        if not is_r1_invalid:
                            misc_output.write(r1)
                        if not is_r2_invalid:
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

    def build_adj(self, matrix, output_prefix):
        self.read1.reset()
        self.read2.reset()
        itr1 = self.read1.fetch(until_eof=True)
        itr2 = self.read2.fetch(until_eof=True)
        print >> sys.stderr, 'populating Hi-C contact matrix with read pairs'
        start_time = time.time()
        total = 0
        for r1 in itr1:
            r2 = itr2.next()
            total += 1
            if total % 5000000 == 0:
                print >> sys.stderr, 'processed %s millions Hi-C contacts.' % str(total/1000000)
            matrix.populate(r1, r2)
        output_filename = output_prefix + ".adj"
        matrix.write(output_filename)
        end_time = time.time()
        print >> sys.stderr, 'cost %s minutes to generate Hi-C contact matrix.' % str(
            round((end_time - start_time) / 60.0, 2))


class SingleRead():
    def __init__(self, read_filename, type='rb'):
        self.read = pysam.AlignmentFile(read_filename, type)
        self.read_filename = read_filename

    def build_bed(self, vector, output_prefix):
        self.read.reset()
        itr = self.read.fetch(until_eof=True)
        print >> sys.stderr, 'populating control bias vector'
        start_time = time.time()
        for r in itr:
            vector.populate(r)
        output_filename = output_prefix + ".bed"
        vector.write(output_filename)
        end_time = time.time()
        print >> sys.stderr, 'cost %s minutes to generate control bias vector.' % str(
            round((end_time - start_time) / 60.0, 2))


def is_unmapped_or_low_mapq(read, mapq):
    return read.is_unmapped or read.mapping_quality < mapq


def is_junk(read1, read2, mapq):
    return is_unmapped_or_low_mapq(read1, mapq) and is_unmapped_or_low_mapq(read2, mapq)


def strand(read):
    """
    Returns:
    reverse strand: +1
    forward strand: -1
    """
    return read.is_reverse * 2 - 1


def is_hic(read1, read2, cutoff, is_r1_invalid, is_r2_invalid):
    if is_r1_invalid or is_r2_invalid:
        return False
    if read1.reference_id != read2.reference_id:
        return True
    else:
        if read1.pos * strand(read1) + read2.pos * strand(read2) > cutoff:
            return True
        else:
            return False


def is_ctl(read1, read2, site, junction, is_read1_invalid, is_read2_invalid):
    is_read1_match = is_match(read1, site)
    is_read2_match = is_match(read2, site)
    if is_read1_match or is_read2_match:
        if is_junction(read1, site, junction) or is_junction(read2, site, junction):
            return False
        elif not is_read1_invalid and not is_read2_invalid:
            if read1.pos * strand(read1) + read2.pos * strand(read2) >= 0:
                return True
            else:
                return False
        else:
            return True
    else:
        return False


def is_match(read, site):
    return regex.match("^%s" % site, read.query_sequence) is not None


def is_junction(read, site, junction):
    if read.cigartuples is None or len(read.cigartuples) == 1:
        return False
    else:
        return junction in read.query_sequence


def is_rlg(read1, read2, is_read1_invalid, is_read2_invalid):
    if is_read1_invalid == is_read2_invalid:
        if read1.pos * strand(read1) + read2.pos * strand(read2) >= 0:
            return True
        else:
            return False
    else:
        return False

