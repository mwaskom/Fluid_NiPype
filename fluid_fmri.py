import os
import re
import sys
import shutil
import argparse
from copy import deepcopy
from datetime import datetime

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.fsl as fsl
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util

from fluid_preproc import preproc
from fluid_source import infosource, datasource

import fluid_utility as flutil

# Parse command line arguments
parser = argparse.ArgumentParser(description="Main interface for GFluid fMRI NiPype code.")

parser.add_argument("-paradigm", metavar="paradigm", 
                    help="experimental paradigm")
parser.add_argument("-subjects", nargs="*",
                    metavar="subject_id",
                    help="run pypeline for subject(s)")
parser.add_argument("-norun",dest="run",action="store_false",
                    help="do not run the pypeline")
parser.add_argument("-inseries",action="store_true",
                    help="force running in series")
args = parser.parse_args()

# Dynamically import the experiment file
# Get the paradigm from the command line, or use the 
# IQ paradigm by default (so this can be importable)
if args.paradigm is not None:
    exp = __import__("%s_experiment" % args.paradigm)
else:
    import iq_experiment as exp
    args.paradigm = "iq"

# Determine the subjects list
# Hierarchy:
# - Command line
# - Subject list defined in experiment module
# - Default subject list, with exclusions from experiment module
if hasattr(exp, "subject_list"):
    subject_list = exp.subject_list
else:
    subject_list = ["gf%02d"%id for id in [4, 5, 9, 13, 14]]
if hasattr(exp, "exclude_subjects"):
    subject_list = [s for s in subject_list if s not in exp.exclude_subjects]

if args.subjects is not None:
    subject_list = args.subjects

# Define some paths
project_dir = "/mindhive/gablab/fluid"
data_dir = os.path.join(project_dir, "Data")
analysis_dir = os.path.join(project_dir, "Analysis/NiPype", args.paradigm)
working_dir = os.path.join(project_dir, "Analysis/NiPype/workingdir", args.paradigm)
report_dir = os.path.join(project_dir, "Analysis/NiPype", args.paradigm, "report")

# Subject source node
# -------------------
subjectsource = pe.Node(util.IdentityInterface(fields=["subject_id"]),
                        iterables = ("subject_id", subject_list),
                        name = "subjectsource")

# Preprocessing
# -------------

# Get preproc input and output
preproc_input = preproc.get_node("inputspec")
preproc_output = preproc.get_node("outputspec")
preproc_report = preproc.get_node("report")

# Preproc datasource node
preprocsource = pe.Node(nio.DataGrabber(infields=["subject_id"],
                                        outfields=["timeseries"],
                                        base_directory=data_dir,
                                        template=exp.source_template,
                                        ),
                        name="preprocsource")
preprocsource.inputs.template_args = exp.template_args

# Preproc Datasink nodes
preprocsink = pe.Node(nio.DataSink(base_directory=analysis_dir),
                      name="preprocsink")

preprocreport = pe.Node(nio.DataSink(base_directory=report_dir),
                        name="preprocreport")

# Preproc filename substitutions
preprocsinksub = pe.Node(util.Merge(len(preproc_output.outputs.__dict__.keys())),
                         name="preprocsinksub")

flutil.get_substitutions(preproc, preproc_output, preprocsinksub)

preprocreportsub = pe.Node(util.Merge(len(preproc_report.outputs.__dict__.keys())),
                              name="preprocreportsub")

flutil.get_substitutions(preproc, preproc_report, preprocreportsub)

# Preproc node substitutions
# NOTE: It would be nice if this were more intuitive, but I haven't
# figured out a good way.  Have to hardcode the node names for now.
sinknodesubs = []
# Shouldn't hurt anything if we just get the maximum number
# of substitutions, although it's a bit messy.
for r in range(4):
    for node in ["art", "dilatemask", "highpass", "masksmoothfunc", "meanfunc2", "realign"]:
        sinknodesubs.append(("_%s%d"%(node, r), "run_%d"%(r+1)))

reportnodesubs = []
for r in range(4):
    for plot in ["displacement", "rotation", "translation"]:
        reportnodesubs.append(("_plot%s%d"%(plot, r), "run_%d"%(r+1)))

# Preproc connections
preproc.connect([
    (subjectsource,    preprocsource, 
        [("subject_id", "subject_id")]),
    (preprocsource,    preproc_input, 
        [("timeseries", "timeseries")]),
    (preproc_output,   preprocsink,   
        [("unsmoothed_timeseries", "preproc.@unsmoothed_timeseries")]),
    (preproc_output,   preprocsink,
        [("smoothed_timeseries", "preproc.@smoothed_timeseries")]),
    (preproc_output,   preprocsink,
        [("example_func", "preproc.@example_func")]),
    (preproc_output,   preprocsink,
        [("functional_mask", "preproc.@functional_mask")]),
    (preproc_output,   preprocsink,
        [("realignment_parameters", "preproc.@realignment_parameters")]),
    (preproc_output,   preprocsink,
        [("outlier_files","preproc.@outlier_files")]),
    (preproc_report,   preprocreport,
        [("rotation_plot", "preproc.@rotation_plot")]),
    (preproc_report,   preprocreport,
        [("translation_plot", "preproc.@translation_plot")]),
    (preproc_report,   preprocreport,
        [("displacement_plot", "preproc.@displacement_plot")]),
    (preprocsinksub,   preprocsink,
        [(("out", lambda x: sinknodesubs + x), "substitutions")]),
    (preprocreportsub, preprocreport,
        [(("out", lambda x: reportnodesubs + x), "substitutions")]),
    ])

# Set up the subject containers
flutil.subject_container(preproc, subjectsource, preprocsink)
flutil.subject_container(preproc, subjectsource, preprocreport)

# Set up the working output
preproc.base_dir = working_dir

# Set the other preprocessing inputs
preproc_input.inputs.hpf_cutoff = exp.hpcutoff/exp.TR
preproc_input.inputs.smooth_fwhm = 5

if __name__ == "__main__":
    if args.run:
        preproc.run(inseries=args.inseries)
