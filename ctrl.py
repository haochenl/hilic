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

    def __init__(self, bam_filename, type='rb'):
        self.bam_reader = pysam.AlignmentFile(bam_filename, type)
        self.bam_filename = bam_filename

    def hic_separate(self, cutoff, mapq, enzyme, output_prefix, threads):
        # the hic reads alignment output
        hic_filename = output_prefix + "_hic.bam"
        hic_output = pysam.AlignmentFile(hic_filename, 'wb', template=self.bam_reader)
        # the control reads alignment output
        ctl_filename = output_prefix + "_hctl.bam"
        ctl_output = pysam.AlignmentFile(ctl_filename, 'wb', template=self.bam_reader)
        # the religation reads alignment output
        rlg_filename = output_prefix + "_rlg.bam"
        rlg_output = pysam.AlignmentFile(rlg_filename, 'wb', template=self.bam_reader)
        # the junk alignment output
        junk_filename = output_prefix + "_hjk.bam"
        junk_output = pysam.AlignmentFile(junk_filename, 'wb', template=self.bam_reader)
        total = 0
        hic_count = 0
        ctl_count = 0
        rlg_count = 0
        junk_count = 0
        print >> sys.stderr, '[process Hi-C alignment files to retrieve contact reads and control reads]'
        start_time = time.time()
        reader = self.bam_reader.fetch(until_eof=True)
        last = None
        for current in reader:
            if last is None:
                last = current
                continue
            if current.query_name == last.query_name:
                total += 1
                is_last_invalid = is_unmapped_or_low_mapq(last, mapq)
                is_current_invalid = is_unmapped_or_low_mapq(current, mapq)
                if is_current_invalid and is_last_invalid:
                    junk_count += 1
                    junk_output.write(last)
                    junk_output.write(current)
                # write the Hi-C output
                elif is_hic(current, last, cutoff, is_current_invalid, is_last_invalid):
                    hic_count += 1
                    hic_output.write(last)
                    hic_output.write(current)
                else:
                    enzyme_site = self._siteDictionary[enzyme]
                    ligation_junction = self._junctionDictionary[enzyme]
                    # write the normalization control output
                    if is_ctl(current, last, enzyme_site, ligation_junction, is_current_invalid, is_last_invalid):
                        ctl_count += 1
                        if not is_last_invalid:
                            ctl_output.write(last)
                        if not is_current_invalid:
                            ctl_output.write(current)
                    # write the re-ligation reads into output
                    elif is_rlg(current, last, is_current_invalid, is_last_invalid):
                        rlg_count += 1
                        rlg_output.write(last)
                        rlg_output.write(current)
                    else:
                        junk_count += 1
                        junk_output.write(last)
                        junk_output.write(current)
                last = None
            else:
                print >> sys.stderr, 'unmatched headers between two reads: (%s, %s)' % (current.query_name, last.query_name)
                last = current
        print >> sys.stderr, 'total number of read pairs processed: %s' % str(total)
        print >> sys.stderr, 'number of Hi-C contacts identified: %s' % str(hic_count)
        print >> sys.stderr, 'number of control read pairs identified: %s' % str(ctl_count)
        print >> sys.stderr, 'number of re-ligation read pairs identified: %s' % str(rlg_count)
        print >> sys.stderr, 'number of junk read pairs: %s' % str(junk_count)
        hic_output.close()
        ctl_output.close()
        rlg_output.close()
        junk_output.close()
        end_time = time.time()
        print >> sys.stderr, 'cost %s minutes to process Hi-C bam files.' % str(
            round((end_time - start_time) / 60.0, 2))
        print >> sys.stderr, '[sort the Hi-C output bam files]'
        start_time = time.time()
        pysam.sort("-@", str(threads), hic_filename, "-o", output_prefix + "_hic_sorted.bam")
        pysam.sort("-@", str(threads), ctl_filename, "-o", output_prefix + "_hctl_sorted.bam")
        pysam.sort("-@", str(threads), rlg_filename, "-o", output_prefix + "_rlg_sorted.bam")
        end_time = time.time()
        print >> sys.stderr, 'cost %s minutes to sort the output bam files.' % str(
            round((end_time - start_time) / 60.0, 2))

    def control_separate(self, cutoff, mapq, enzyme, output_prefix, threads):
        ctl_filename = output_prefix + "_xctl.bam"
        ctl_output = pysam.AlignmentFile(ctl_filename, 'wb', template=self.bam_reader)
        misc_filename = output_prefix + "_misc.bam"
        misc_output = pysam.AlignmentFile(misc_filename, 'wb', template=self.bam_reader)
        junk_filename = output_prefix + "_xjk.bam"
        junk_output = pysam.AlignmentFile(junk_filename, 'wb', template=self.bam_reader)
        total = 0
        ctl_count = 0
        misc_count = 0
        junk_count = 0
        print >> sys.stderr, '[process control alignment files to retrieve control reads]'
        start_time = time.time()
        reader = self.bam_reader.fetch(until_eof=True)
        last = None
        for current in reader:
            if last is None:
                last = current
                continue
            if last.query_name == current.query_name:
                total += 1
                is_last_invalid = is_unmapped_or_low_mapq(last, mapq)
                is_current_invalid = is_unmapped_or_low_mapq(current, mapq)
                # write into junk output if read pair meets Hi-C conditions
                if is_hic(last, current, cutoff, is_last_invalid, is_current_invalid):
                    junk_count += 1
                    junk_output.write(last)
                    junk_output.write(current)
                else:
                    enzyme_site = self._siteDictionary[enzyme]
                    ligation_junction = self._junctionDictionary[enzyme]
                    # write the normalization control output
                    if is_ctl(last, current, enzyme_site, ligation_junction, is_last_invalid, is_current_invalid):
                        ctl_count += 1
                        if not is_last_invalid:
                            ctl_output.write(last)
                        if not is_current_invalid:
                            ctl_output.write(current)
                    # write the miscellaneous control read to output
                    else:
                        misc_count += 1
                        if not is_last_invalid:
                            misc_output.write(last)
                        if not is_current_invalid:
                            misc_output.write(current)
                last = None
            else:
                print >> sys.stderr, 'unmatched headers between two reads: (%s, %s)' % (last.query_name, current.query_name)
                last = current
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
        print >> sys.stderr, '[sort the control output bam files]'
        start_time = time.time()
        pysam.sort("-@", str(threads), ctl_filename, "-o", output_prefix + "_xctl_sorted.bam")
        pysam.sort("-@", str(threads), misc_filename, "-o", output_prefix + "_misc_sorted.bam")
        end_time = time.time()
        print >> sys.stderr, 'cost %s minutes to sort the output bam files.' % str(
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
        if abs(read1.pos - read2.pos) > cutoff:
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

