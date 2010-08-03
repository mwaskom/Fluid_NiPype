"""
    Nipype interface script for the Fluid Intelligence project.
"""


import os
import sys
from datetime import datetime

import nipype.pipeline.engine as pe

from fluid_preproc import preproc
from fluid_datasource import infosource, datasource

if not pe.__file__.startswith("/software/python/nipype0.3"):
    sys.exit("ERROR: Not using nipype0.3")

subject_list = ["SMARTER_SP14"]

par = "nback"

highpass = preproc.get_node("highpass")
if par == "nback":
    info = dict(func=[["subj_id", "nii", "NBack1"]],
                struct=[["subj_id", "nii", "mprage"]])
    hpcutoff = 120
    TR = 3.

highpass.inputs.suffix = "_hpf"
highpass.inputs.op_string = "-bptf %s -1" % (hpcutoff/TR)

smooth = preproc.get_node("smooth")
smooth.inputs.fwhm = 6.

datasource.inputs.template_args = info

infosource.iterables = ("subj_id", subject_list)

inputnode = preproc.get_node("inputspec")

preproc.connect([(infosource, datasource, [("subj_id", "subj_id")]),
                 (datasource, inputnode, [("struct", "struct"),
                                          ("func", "func")])])

datestamp = str(datetime.now())[:10]
codepath = os.path.split(os.path.abspath(__file__))[0]
crashdir = os.path.abspath("%s/crashdumps/%s" % (codepath, datestamp))

if not os.path.isdir(crashdir):    
    os.makedirs(crashdir)

preproc.config = dict(crashdump_dir=crashdir) 
preproc.base_dir = os.path.abspath("../nipype_output/workingdir/")

if __name__ == "__main__":
    preproc.run()
    preproc.write_graph(graph2use = "flat")
