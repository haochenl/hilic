#!/usr/bin/env python

"""
The hilic file I/O utilities
"""

import sys
import os

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
