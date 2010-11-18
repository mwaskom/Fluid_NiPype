import os
import re
import sys
import imp
import shutil
import argparse
import inspect
from copy import deepcopy
from datetime import datetime

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.fsl as fsl
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util
from nipype.interfaces.base import Bunch

from workflows.preproc import preproc
from workflows.registration import registration
from workflows.fsl_model import fsl_model
#from workflows.fluid_fixed_fx import fixed_fx

import fluid_utility as flutil

# Parse command line arguments
parser = argparse.ArgumentParser(description="Main interface for GFluid fMRI NiPype code.")

parser.add_argument("-paradigm", metavar="paradigm", 
                    help="experimental paradigm")
parser.add_argument("-subjects", nargs="*",
                    metavar="subject_id",
                    help="process subject(s)")
parser.add_argument("-workflows", nargs="*",
                    metavar="WF",
                    help="which workflows to run")
parser.add_argument("-norun",dest="run",action="store_false",
                    help="do not run the pypeline")
parser.add_argument("-inseries",action="store_true",
                    help="force running in series")
args = parser.parse_args()

# Define a default paradigm, for importability and convenince
default_paradigm = "iq"

# Dynamically import the experiment file
# Get the paradigm from the command line, or use the 
# IQ paradigm by default (so this can be importable)
if args.paradigm is None:
    args.paradigm = default_paradigm
try:   
    exp = __import__("experiments.%s_experiment"%args.paradigm,
                     fromlist=["experiments"])
except ImportError:
    exp = __import__("%s_experiment"%args.paradigm)

if args.workflows is None:
    args.workflows = []
elif args.workflows == ["all"]:
    args.workflows = ["preproc", "reg", "model"]

# Determine the subjects list
# Hierarchy:
# - Command line
# - Subject list defined in experiment module
# - Default subject list, with exclusions from experiment module
if args.subjects is not None:
    subject_list = args.subjects
elif hasattr(exp, "subject_list"):
    subject_list = exp.subject_list
else:
    subject_list = ["gf%02d"%id for id in [4, 5, 9, 13, 14]]

if hasattr(exp, "exclude_subjects"):
    subject_list = [s for s in subject_list if s not in exp.exclude_subjects]


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
for r in range(exp.nruns):
    for node in ["art", "dilatemask", "highpass", "masksmoothfunc", 
                 "extractref", "meanfunc2", "realign"]:
        preprocsinknodesubs.append(("_%s%d"%(node, r), "run_%d"%(r+1)))

preprocreportnodesubs = []
for r in range(exp.nruns):
    for plot in ["displacement", "rotation", "translation"]:
        preprocreportnodesubs.append(("_plot%s%d"%(plot, r), "run_%d"%(r+1)))
    for img in ["example", "mean"]:
        preprocreportnodesubs.append(("_%sslice%d"%(img, r), "run_%d"%(r+1)))

# Preproc connections
preproc.connect([
    (subjectsource,    preprocsource, 
        [("subject_id", "subject_id")]),
    (preprocsinksub,   preprocsink,
        [(("out", lambda x: preprocsinknodesubs + x), "substitutions")]),
    (preprocreportsub, preprocreport,
        [(("out", lambda x: preprocreportnodesubs + x), "substitutions")]),
    ])

# Input connections
flutil.connect_inputs(preproc, preprocsource, preproc_input)

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
                                    outfields=["warpfield", "mean_func",
                                               "vol_timeseries", "surf_timeseries"],
                                    base_directory = project_dir,
                                    sort_filelist=True),
                    name="regsource")

regsource.inputs.template = "Analysis/NiPype/" + args.paradigm + "/%s/preproc/run_?/%s.nii.gz"
regsource.inputs.field_template = dict(warpfield="Data/%s/registration/%s.nii.gz")
regsource.inputs.template_args = dict(mean_func=[["subject_id", "mean_func"]],
                                      example_func=[["subject_id", "example_func"]],
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
for r in range(exp.nruns):
    for node in ["func2anat","warptimeseries"]:
        regsinknodesubs.append(("_%s%d"%(node, r), "run_%d"%(r+1)))

regreportnodesubs = []
for r in range(exp.nruns):
    for node in ["exfuncwarppng", "func2anat", "func2anatpng"]:
        regreportnodesubs.append(("_%s%d"%(node, r), "run_%d"%(r+1)))

# Registration connections
registration.connect([
    (subjectsource, regsource,  
        [("subject_id", "subject_id")]),
    (subjectsource, reg_input,
        [("subject_id", "subject_id")]),
    (regsinksub,   regsink,
        [(("out", lambda x: regsinknodesubs + x), "substitutions")]),
    (regreportsub, regreport,
        [(("out", lambda x: regreportnodesubs + x), "substitutions")]),
    ])

# Connect the inputs
flutil.connect_inputs(registration, regsource, reg_input)

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

# First-Level Model
# -----------------

# Moderate hack to strip deepcopy-offensive stuff from an experiment module
# Fixes bug in Python(!)
class Foo(object): pass
experiment = Foo()
for k,v in exp.__dict__.items():
    if not k.startswith("__") and not inspect.ismodule(v):
        setattr(experiment, k, v)

# Model relevant functions
def count_runs(runs):
    """Count the number of functional timeseries globbed for each subject."""
    if not isinstance(runs, list):
        runs = [runs]
    return len(runs)

def get_contrast_idx(contrast_name):
    """Return the index corresponding to the name of a contrast."""
    return [i for i, data in enumerate(exp.contrasts) if data[0] == contrast_name]

def subjectinfo(info, exp):
    """Return a subjectinfo list of bunches.  

    Parameters
    ----------
    info : list
        [subject_id, nruns]
    exp : experiment object
        Basically an experiment module stipped of deepcopy-offensive stuff

    """
    subject_id = info[0]
    nruns = info[1]
    output = []
    events = exp.events
    for r in range(nruns):
        run_number = r+1
        onsets = []
        durations = []        
        for event in events:

            # Get a path to the parfile for this event
            vars = locals()
            arg_tuple = tuple([vars[arg] for arg in exp.parfile_args])
            parfile = os.path.join(exp.parfile_base_dir, exp.parfile_template%arg_tuple)
            
            # Get the event onsets and durations from the parfile
            o, d = flutil.parse_par_file(parfile)
            onsets.append(o)
            durations.append(d)

        output.insert(r,
                      Bunch(conditions=events,
                            onsets=deepcopy(onsets),
                            durations=deepcopy(durations),
                            amplitudes=None,
                            tmod=None,
                            pmod=None,
                            regressor_names=None,
                            regressors=None))
    return output

# Set the template and args for the timeseries in each space
timeseries_template = dict(volume = dict(
    timeseries=os.path.join("Analysis/NiPype",args.paradigm,"%s/registration/run_?/%s.nii.gz")))

timeseries_args = dict(volume= [["subject_id", "warped_timeseries"]])

# Set up the model workflow for each space 

model = fsl_model.clone("volume_model")

# Get model input and output
model_input = model.get_node("inputspec")
model_output = model.get_node("outputspec")
model_report = model.get_node("report")

# Set up the date sources for each space
modelsource = pe.Node(nio.DataGrabber(infields=["subject_id"],
                                      outfields=["timeseries",
                                                 "outlier_files",
                                                 "realignment_parameters"],
                                      base_directory = project_dir,
                                      sort_filelist=True),
                      name="modelsource")

modelsource.inputs.template = os.path.join(
    "Analysis/NiPype",args.paradigm ,"%s/preproc/run_?/%s")

modelsource.inputs.template_args = dict(
    outlier_files=[["subject_id", "outlier_files.txt"]],
    realignment_parameters=[["subject_id", "realignment_parameters.par"]])
modelsource.inputs.field_template = timeseries_template["volume"]
modelsource.inputs.template_args["timeseries"] = timeseries_args["volume"]

# Model Datasink nodes
modelsink = pe.Node(nio.DataSink(base_directory=analysis_dir),
                  name="modelsink")

modelreport = pe.Node(nio.DataSink(base_directory=report_dir),
                    name="modelreport")

modelreportsub = pe.Node(util.Merge(len(model_report.outputs.__dict__.keys())),
                       name="modelreportsub")

flutil.get_substitutions(model, model_report, modelreportsub)

# Model node substitutions
# NOTE: It would be nice if this were more intuitive, but I haven't
# figured out a good way.  Have to hardcode the node names for now.
modelsinknodesubs = []
# Shouldn't hurt anything if we just get the maximum number
# of substitutions, although it's a bit messy.
for r in range(exp.nruns):
    for node in ["modelestimate"]:
        modelsinknodesubs.append(("_%s%d"%(node, r), "run_%d"%(r+1)))
modelsink.inputs.substitutions = modelsinknodesubs

modelreportnodesubs = [("_contrast_","stats/")]
for r in range(exp.nruns):
    for node in ["modelgen", "sliceresidual", "slicestats"]:
        modelreportnodesubs.append(("_%s%d"%(node, r), "run_%d"%(r+1)))

# Define a node to seed the contrasts iterables with the contrast name
connames = pe.Node(util.IdentityInterface(fields=["contrast"]),
                   iterables = ("contrast", [c[0] for c in exp.contrasts]),
                   name="contrastnames")

# Get a reference to the selectcontrast node
selectcontrast = model.get_node("selectcontrast")

# Define a node to count the number of functional runs
# and merge that with the subject_id to pass to the 
# subjectinfo function
runinfo = pe.Node(util.Merge(2), name="runinfo")

# Model connections
model.connect([
    (subjectsource, modelsource,  
        [("subject_id", "subject_id")]),
    (subjectsource, runinfo,
        [("subject_id", "in1")]),
    (modelsource, runinfo,
        [(("timeseries", count_runs), "in2")]),
    (runinfo, model_input,
        [(("out", subjectinfo, experiment), "subject_info")]),
    (connames, selectcontrast,
        [(("contrast", get_contrast_idx), "index")]),
    (modelreportsub, modelreport,
        [(("out", lambda x: modelreportnodesubs + x), "substitutions")]),
    ])

# Connect inputs
flutil.connect_inputs(model, modelsource, model_input)

# Set up the subject containers
flutil.subject_container(model, subjectsource, modelsink)
flutil.subject_container(model, subjectsource, modelreport)

# Connect the heuristic outputs to the datainks
flutil.sink_outputs(model, model_output, modelsink, "model.%s"%"volume")
flutil.sink_outputs(model, model_report, modelreport, "model.%s"%"volume")

# Set the other model inputs
model_input.inputs.TR = exp.TR
model_input.inputs.units = exp.units
model_input.inputs.hpf_cutoff = exp.hpcutoff
model_input.inputs.HRF_bases = exp.fsl_bases
model_input.inputs.contrasts = exp.contrasts

# Volume model working output
model.base_dir = working_dir

# Set crashdumps
flutil.archive_crashdumps(model)

# Across-run Fixed Fx



if __name__ == "__main__":
    if args.run:
        if __file__ == "fluid_fmri.py":
            report_script = "python /mindhive/gablab/fluid/NiPype_Code/reporting/build_report.py "
        
        if "preproc" in args.workflows:
            preproc.run(inseries=args.inseries)
            os.system(report_script + " ".join(subject_list))
        
        if "reg" in args.workflows:
            registration.run(inseries=args.inseries)
            os.system(report_script + " ".join(subject_list))
        
        if "model" in args.workflows:
            model.run(inseries=args.inseries)
            os.system(report_script + " ".join(subject_list))
