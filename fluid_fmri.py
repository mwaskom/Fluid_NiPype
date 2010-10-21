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
from fluid_registration import registration

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
                                        sort_filelist=True,
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
preprocsinknodesubs = []
# Shouldn't hurt anything if we just get the maximum number
# of substitutions, although it's a bit messy.
for r in range(4):
    for node in ["art", "dilatemask", "highpass", "masksmoothfunc", "meanfunc2", "realign"]:
        preprocsinknodesubs.append(("_%s%d"%(node, r), "run_%d"%(r+1)))

preprocreportnodesubs = []
for r in range(4):
    for plot in ["displacement", "rotation", "translation"]:
        preprocreportnodesubs.append(("_plot%s%d"%(plot, r), "run_%d"%(r+1)))

# Preproc connections
preproc.connect([
    (subjectsource,    preprocsource, 
        [("subject_id", "subject_id")]),
    (preprocsource,    preproc_input, 
        [("timeseries", "timeseries")]),
    (preprocsinksub,   preprocsink,
        [(("out", lambda x: preprocsinknodesubs + x), "substitutions")]),
    (preprocreportsub, preprocreport,
        [(("out", lambda x: preprocreportnodesubs + x), "substitutions")]),
    ])
# Set up the subject containers
flutil.subject_container(preproc, subjectsource, preprocsink)
flutil.subject_container(preproc, subjectsource, preprocreport)

# Connect the heuristic outputs to the datainks
flutil.sink_outputs(preproc, preproc_output, preprocsink, "preproc")
flutil.sink_outputs(preproc, preproc_report, preprocreport, "preproc")

# Set up the working output
preproc.base_dir = working_dir

# Archive crashdumps
flutil.archive_crashdumps(preproc)

# Set the other preprocessing inputs
preproc_input.inputs.hpf_cutoff = exp.hpcutoff/exp.TR
preproc_input.inputs.smooth_fwhm = 5


# Registration
# ------------

# Get registration input and output
reg_input = registration.get_node("inputspec")
reg_output = registration.get_node("outputspec")
reg_report = registration.get_node("report")

# Registration datasource node
regsource = pe.Node(nio.DataGrabber(infields=["subject_id"],
                                    outfields=["example_func", "warpfield",
                                               "vol_timeseries", "surf_timeseries"],
                                    base_directory = project_dir,
                                    sort_filelist=True),
                    name="regsource")

regsource.inputs.template = "Analysis/NiPype/" + args.paradigm + "/%s/preproc/run_?/%s.nii.gz"
regsource.inputs.field_template = dict(warpfield="Data/%s/registration/%s.nii.gz")
regsource.inputs.template_args = dict(example_func=[["subject_id", "example_func"]],
                                      vol_timeseries=[["subject_id", "smoothed_timeseries"]],
                                      surf_timeseries=[["subject_id", "unsmoothed_timeseries"]],
                                      warpfield=[["subject_id","warpfield"]])

# Registration Datasink nodes
regsink = pe.Node(nio.DataSink(base_directory=analysis_dir),
                  name="regsink")

regreport = pe.Node(nio.DataSink(base_directory=report_dir),
                    name="regreport")

# Registration filename substitutions
regsinksub = pe.Node(util.Merge(len(reg_output.outputs.__dict__.keys())),
                     name="regsinksub")

flutil.get_substitutions(registration, reg_output, regsinksub)

regreportsub = pe.Node(util.Merge(len(reg_report.outputs.__dict__.keys())),
                       name="regreportsub")

flutil.get_substitutions(registration, reg_report, regreportsub)

# Registration node substitutions
# NOTE: It would be nice if this were more intuitive, but I haven't
# figured out a good way.  Have to hardcode the node names for now.
regsinknodesubs = []
# Shouldn't hurt anything if we just get the maximum number
# of substitutions, although it's a bit messy.
for r in range(4):
    for node in ["func2anat","warptimeseries"]:
        regsinknodesubs.append(("_%s%d"%(node, r), "run_%d"%(r+1)))

regreportnodesubs = []
for r in range(4):
    for node in ["exfuncwarppng", "func2anat", "func2anatpng"]:
        regreportnodesubs.append(("_%s%d"%(node, r), "run_%d"%(r+1)))

# Registration connections
registration.connect([
    (subjectsource, regsource,  
        [("subject_id", "subject_id")]),
    (subjectsource, reg_input,
        [("subject_id", "subject_id")]),
    (regsource,     reg_input,  
        [("vol_timeseries", "vol_timeseries"),
         ("surf_timeseries", "surf_timeseries"),
         ("example_func", "example_func"),
         ("warpfield", "warpfield")]),
    (regsinksub,   regsink,
        [(("out", lambda x: regsinknodesubs + x), "substitutions")]),
    (regreportsub, regreport,
        [(("out", lambda x: regreportnodesubs + x), "substitutions")]),
    ])

# Set up the subject containers
flutil.subject_container(registration, subjectsource, regsink)
flutil.subject_container(registration, subjectsource, regreport)

# Connect the heuristic outputs to the datainks
flutil.sink_outputs(registration, reg_output, regsink, "registration")
flutil.sink_outputs(registration, reg_report, regreport, "registration")

# Registration working output
registration.base_dir = working_dir

# Set crashdumps
flutil.archive_crashdumps(registration)

if __name__ == "__main__":
    if args.run:
        
        preproc.run(inseries=args.inseries)
        for subj in subject_list:
            os.system("python /mindhive/gablab/fluid/NiPype_Code/build_subj_report.py %s"%subj)
        os.system("python /mindhive/gablab/fluid/NiPype_Code/build_report_homepage.py")
        
        registration.run(inseries=args.inseries)
        for subj in subject_list:
            os.system("python /mindhive/gablab/fluid/NiPype_Code/build_subj_report.py %s"%subj)
        os.system("python /mindhive/gablab/fluid/NiPype_Code/build_report_homepage.py")
