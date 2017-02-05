#!/usr/bin/env python

"""
The hilic threading utilities
"""

import subprocess

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
