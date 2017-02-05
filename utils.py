#!/usr/bin/env python

"""
The hilic common utilities
"""

import sys
import os


def count_complete_process(processes):
    """
    Moniter a list of processes and return the number of running processes
    :param processes: a list of processes
    :return: the number of running processes
    """
    count = 0
    for proc, filename in processes:
        status = proc.poll()
        if status is None:
            continue
        else:
            count += 1
    return count


def check_file_status(filename_list):
    """
    Check a list of file names and exit if not a regular file
    :param filename_list: list of file names
    :return:
    """
    for filename in filename_list:
        if os.path.isdir(filename):
            print >> sys.stderr, 'file name is a directory not a regular file: %s' % str(filename)
            print >> sys.stderr, 'directory path not supported in the config file'
            sys.exit(1)
        else:
            if not os.path.exists(filename):
                print >> sys.stderr, 'listed file does not exist: %s' % str(filename)
                sys.exit(1)
