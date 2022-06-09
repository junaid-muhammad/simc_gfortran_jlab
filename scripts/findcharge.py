#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2022-06-04 13:11:38 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import re
import sys, os

################################################################################################################################################
'''
User Inputs
'''
ROOTPrefix = sys.argv[1]
runNum = sys.argv[2]
MaxEvent=sys.argv[3]

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
import ltsep as lt 

# Add this to all files for more dynamic pathing
USER =  lt.SetPath(os.path.realpath(__file__)).getPath("USER") # Grab user info for file finding
HOST = lt.SetPath(os.path.realpath(__file__)).getPath("HOST")
REPLAYPATH = lt.SetPath(os.path.realpath(__file__)).getPath("REPLAYPATH")
UTILPATH = lt.SetPath(os.path.realpath(__file__)).getPath("UTILPATH")
ANATYPE=lt.SetPath(os.path.realpath(__file__)).getPath("ANATYPE")

################################################################################################################################################

# Open report file to grab prescale values and tracking efficiency
report = UTILPATH+"/REPORT_OUTPUT/Analysis/General/%s_%s_%s.report" % (ROOTPrefix,runNum,MaxEvent)

with open(report) as f:
    # Search for keywords, then save as value in dictionary
    for line in f:
        data = line.split(':')
        if 'BCM1_Beam_Cut_Charge' in data[0]:
            charge = int(re.sub("\D","","%s" % data[1]))
print(charge)