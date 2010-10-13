# vi: set ft=python sts=4 ts=4 sw=4 et:
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
"""
    Nipype interface script for the Fluid Intelligence project.
    EDITS for rush analysis pre DARPA meeting (volume only)
"""

import os
import re
import sys
import shutil
from copy import deepcopy
from datetime import datetime

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.spm as spm
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util

# Catch when we"re not using the right environment
if (not pe.__file__.startswith("/software/python/nipype0.3") 
    and not pe.__file__.startswith("/u2/mwaskom/nipype")):
    sys.exit("ERROR: Not using nipype0.3")

import argparse

from new_preproc import preproc
from fluid_fsl_model import fsl_vol_model
from fluid_fixed_fx import vol_fixed_fx as fixed_fx

import fluid_utility_funcs as fuf
if __name__ == "__main__":
    """ Handle command line arguments to control the analysis """
    parser = argparse.ArgumentParser(description="Main interface for GFluid NiPype code.")
    parser.add_argument("-paradigm", metavar="paradigm", required=True,
                        help="experimental paradigm")
    parser.add_argument("-subject", dest="subjects", 
                        metavar="subject_id", action="append",
                        help="run pypeline for subject(s)")
    parser.add_argument("-norun",dest="run",action="store_false",
                        help="do not run the pypeline")
    parser.add_argument("-nograph",dest="write_graph",action="store_false",
                        help="do not write a graph")
    parser.add_argument("-inseries",action="store_true",
                        help="force running in series")
    parser.add_argument("-runspm",action="store_true",
                        help="run spm first level")
    args = parser.parse_args()
    """ Dynamically import the experiment file """
else:
    class Foo(): pass
    args = Foo()
    args.paradigm = "nback"
    args.run = False
    args.write_graph = False
exp = __import__("%s_experiment" % args.paradigm)

if args.paradigm == "nback":
	subject_list = ["gf%02d"%id for id in [4,5,9,13,14,18,20,27]]
elif args.paradigm == "iq":
    subject_list = ["gf%02d"%id for id in [4,5,9,13,14,18,20,21,23,27]]
elif args.paradigm == "mot_block":
    subject_list = ["gf%02d"%id for id in [4,18,21,27]]
elif args.paradigm == "mot_jitter":
    subject_list = ["gf%02d"%id for id in [4,18,27]]
if args.subjects:
    subject_list = args.subjects
print "Subjects:"," ".join(subject_list)

firstlevel = pe.Workflow(name="level1")

# Data and subject source
data_dir = "/mindhive/gablab/fluid/Data"
exp.data_dir = data_dir

infosource = pe.Node(interface=util.IdentityInterface(fields=["subject_id"]),
                     name="infosource")

infosource.iterables = ("subject_id", subject_list)

datasource = pe.Node(interface=nio.DataGrabber(infields=["subject_id"],
                                               outfields=["func", "struct"],
                                               base_directory = data_dir),
                     name="datasource")

datasource.inputs.template = exp.source_template
datasource.inputs.template_args = exp.template_args

firstlevel.connect([
    (infosource, datasource, [("subject_id", "subject_id")]),
    (infosource, preproc, [("subject_id", "inputspec.subject_id")]),
    (datasource, preproc, [("func", "inputspec.func")])
    ])

contrasts = exp.contrasts

fsl_vol_model.inputs.modelspec.time_repetition = exp.TR
fsl_vol_model.inputs.modelspec.high_pass_filter_cutoff = exp.hpcutoff
fsl_vol_model.inputs.modelspec.input_units = exp.units
fsl_vol_model.inputs.modelspec.output_units = exp.units

fsl_vol_model.inputs.level1design.contrasts = contrasts
fsl_vol_model.inputs.level1design.interscan_interval = exp.TR
fsl_vol_model.inputs.level1design.bases = exp.fsl_bases

selectvolcontrast = fsl_vol_model.get_node("selectcontrast")
selectvolcontrast.iterables = ('index', [[i] for i in range(len(contrasts))])

def sort_copes(files):
    numelements = len(files[0])
    outfiles = []
    for i in range(numelements):
        outfiles.insert(i,[])
        for j, elements in enumerate(files):
            outfiles[i].append(elements[i])
    return outfiles

# Fuck you, deepcopy()
class Foo(): pass
info = Foo()
for k, v in exp.__dict__.items():
    if not k.startswith("__") and not k == "re":
        setattr(info, k, v)

firstlevel.connect([
    (infosource, fsl_vol_model, 
        [(("subject_id", fuf.subjectinfo, info),"modelspec.subject_info")]),
    (preproc, fsl_vol_model, 
        [("warpfunc.out_file", "modelspec.functional_runs"),
         ("art.outlier_files", "modelspec.outlier_files"),
         ("realign.par_file", "modelspec.realignment_parameters"),
         ("warpfunc.out_file","modelestimate.in_file")]),
      (fsl_vol_model, fixed_fx,
        [(("contrastestimate.copes", sort_copes), "copemerge.in_files"),
         (("contrastestimate.varcopes", sort_copes),"varcopemerge.in_files"),
         (("contrastestimate.copes", lambda x: len(x)), "l2model.num_copes")]),
    ])

# File crashdumps by date
datestamp = str(datetime.now())[:10]
codepath = os.path.split(os.path.abspath(__file__))[0]
crashdir = os.path.abspath("%s/crashdumps/%s" % (codepath, datestamp))
if not os.path.isdir(crashdir):    
    os.makedirs(crashdir)
firstlevel.config = dict(crashdump_dir=crashdir) 

""" Setup the output """
working_output = os.path.join(
    os.path.abspath("../working"), args.paradigm)
firstlevel.base_dir = working_output

output_base = os.path.join(os.path.abspath(".."), args.paradigm)

if args.run:
    firstlevel.run(inseries=args.inseries)
if args.write_graph:
    firstlevel.write_graph(graph2use="flat")


