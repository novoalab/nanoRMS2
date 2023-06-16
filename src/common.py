#!/usr/bin/env python3

import os, sys, resource
from datetime import datetime

VERSION = '2.0a'

def make_counter(count=0):
    """Simple counter"""
    def inner():
        nonlocal count
        count += 1
        return count
    return inner

def warning(info, m=3, counter=make_counter()):
    """Return at most m of warning messages"""
    # update counter
    count = counter()#; print(count)
    # skip if too many warnings returned
    if m and count>m: return
    if info: logger(info, add_timestamp=0, add_memory=0)
    if count==m:
        logger("Future warnings will be skipped!", add_timestamp=0, add_memory=0)

def memory_usage(childrenmem=True, div=1024.):
    """Return memory usage in MB including children processes"""
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / div
    if childrenmem:
        mem += resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / div
    return mem
        
def logger(info, add_timestamp=1, add_memory=1, out=sys.stderr):
    """Report nicely formatted stream to stderr"""
    info = info.rstrip('\n')
    memory = timestamp = ""
    if add_timestamp:
        timestamp = "[%s]"%str(datetime.now()).split(".")[0] #"[%s]"%datetime.ctime(datetime.now())
    if add_memory:
        memory = " [mem: %5.0f MB]"%memory_usage()
    out.write("%s %s%s\n"%(timestamp, info, memory))

