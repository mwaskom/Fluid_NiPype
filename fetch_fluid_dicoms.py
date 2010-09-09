#! /usr/bin/env python
import sys
import os

if len(sys.argv) != 3:
    sys.exit("USAGE: fetch_fluid_dicoms subj type")
else:
    subj = sys.argv[1]
    type = sys.argv[2]

targdir = "/mindhive/gablab/fluid/Data/%s/dicom/%s"%(subj,type)

os.system("fetch_dicoms -s %s -l -d %s"%(subj,targdir))
