#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2022-11-29 12:00:16 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import sys

inp_f = sys.argv[1]

with open(inp_f, 'r') as f:
    f_data = f.read()
    
print(f_data)