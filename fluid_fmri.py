#! /usr/bin/env python
"""
Main interface for GFluid fMRI Nipype code.
"""
import os
import argparse
import inspect
from copy import deepcopy
from os.path import join as pjoin

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
from nipype.interfaces.base import Bunch

from workflows.preproc import get_preproc_workflow
from workflows.fsl_model import get_model_workflow
from workflows.registration import get_registration_workflow
from workflows.fixed_fx import get_fixedfx_workflow

import fluid_utility as flutil

# Parse command line arguments
parser = argparse.ArgumentParser(description="Main interface for GFluid fMRI Nipype code.")

parser.add_argument("-paradigm", metavar="paradigm", 
                    help="experimental paradigm")
parser.add_argument("-subjects", nargs="*",
                    metavar="subject_id",
                    help="process subject(s)")
parser.add_argument("-workflows", nargs="*",
                    metavar="WF",
                    help="which workflows to run")
parser.add_argument("-modelspace", metavar="space", dest="space", default="volume",
                    help="space to run the model in (volume, surface)")
parser.add_argument("-norun",dest="run",action="store_false",
                    help="do not run the pypeline")
parser.add_argument("-inseries",action="store_true",
                    help="force running in series")
args = parser.parse_args()

# Define a default paradigm, for importability and convenince
default_paradigm = "iq"

# Dynamically import the experiment file
if args.paradigm is None:
    args.paradigm = default_paradigm
# Look for paradgim.py in an experiments package
try:   
    exp = __import__("experiments." + args.paradigm,
                     fromlist=["experiments"])
# Or maybe just get it from the current directory
except ImportError:
    exp = __import__(args.paradigm)

if args.workflows is None:
    args.workflows = []
elif args.workflows == ["all"]:
    args.workflows = ["preproc", "reg", "model", "ffx"]

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

print "Subjects: ", " ".join(subject_list)

# Define some paths
project_dir = "/mindhive/gablab/fluid"
data_dir = os.path.join(project_dir, "Data")
analysis_dir = os.path.join(project_dir, "Analysis/Nipype", args.paradigm)
working_dir = os.path.join(project_dir, "Analysis/Nipype/workingdir", args.paradigm)

# Set the Freesurfer Subjects Directory 
os.environ["SUBJECTS_DIR"] = data_dir

# Subject source node
# -------------------
subjectsource = pe.Node(util.IdentityInterface(fields=["subject_id"]),
                        iterables = ("subject_id", subject_list),
                        name = "subjectsource")

# Preprocessing
# =============

# Get preproc input and output
preproc, preproc_input, preproc_output = get_preproc_workflow(name="preproc",
                                                              b0_unwarp=False,
                                                              anat_reg=True)

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

# Preproc filename substitutions
preprocsub = pe.Node(util.Merge(len(preproc_output.outputs.__dict__)),
                         name="preprocsinksub")

flutil.get_output_substitutions(preproc, preproc_output, preprocsub)

# Preproc node substitutions
# NOTE: It would be nice if this were more intuitive, but I haven't
# figured out a good way.  Have to hardcode the node names for now.
preproc_mapnodes = ["art", "dilatemask", "highpass", "masksmoothfunc", "extractref", "tsmovie",
                    "meanfunc3", "realign", "plotmean", "func2anat", "func2anatpng"]
for plot in ["displacement", "rotation", "translation"]:
    preproc_mapnodes.append("plot%s"%plot)
for img in ["example", "mean"]:
    preproc_mapnodes.append("%sslice"%img)
preprocsinksubs = flutil.get_mapnode_substitutions(exp.nruns, preproc_mapnodes)


flutil.set_substitutions(preproc, preprocsink, preprocsub, preprocsinksubs)

# Preproc connections
preproc.connect(subjectsource, "subject_id", preprocsource, "subject_id")
preproc.connect(subjectsource, "subject_id", preproc_input, "subject_id")

# Input connections
flutil.connect_inputs(preproc, preprocsource, preproc_input)

# Set up the subject containers
flutil.subject_container(preproc, subjectsource, preprocsink)

# Connect the heuristic outputs to the datainks
flutil.sink_outputs(preproc, preproc_output, preprocsink, "preproc")

# Set up the working output
preproc.base_dir = working_dir

# Archive crashdumps
flutil.archive_crashdumps(preproc)

# Set the other preprocessing inputs
preproc_input.inputs.hpf_cutoff = exp.hpcutoff/exp.TR
preproc_input.inputs.smooth_fwhm = 5


# First-Level Model
# =================

# Moderate hack to strip deepcopy-offensive stuff from an experiment module
# Fixes bug in Python(!)
class Experiment(object): pass
experiment = Experiment()
for k,v in exp.__dict__.items():
    if not k.startswith("__") and not inspect.ismodule(v):
        setattr(experiment, k, v)

# Model relevant functions
def count_runs(runs):
    """Count the number of functional timeseries being analyzed for each subject."""
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
    if subject_id.endswith("p"):
        day = 2
    else:
        day = 1
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
                            regressor_names=deepcopy(exp.events),
                            amplitudes=None,
                            tmod=None,
                            pmod=None,
                            regressors=None))
    return output



# Map space to smoothing
spacedict = dict(volume="smoothed", surface="unsmoothed")

# Get model input and output
model, model_input, model_output = get_model_workflow(name="%s_model"%spacedict[args.space])

# Set up the date sources for each space
modelsource = pe.Node(nio.DataGrabber(infields=["subject_id"],
                                      outfields=["timeseries",
                                                 "outlier_files",
                                                 "overlay_background",
                                                 "realignment_parameters"],
                                      base_directory = project_dir,
                                      sort_filelist=True),
                      name="modelsource")

modelsource.inputs.template = os.path.join(analysis_dir, "%s/%s/run_?/%s")

modelsource.inputs.template_args = dict(
    outlier_files=[["subject_id", "preproc", "outlier_files.txt"]],
    realignment_parameters=[["subject_id", "preproc", "realignment_parameters.par"]],
    overlay_background=[["subject_id", "preproc", "example_func.nii.gz"]],
    timeseries=[["subject_id", "preproc", "%s_timeseries.nii.gz"%spacedict[args.space]]])


# Model Datasink nodes
modelsink = pe.Node(nio.DataSink(base_directory=analysis_dir),
                  name="modelsink")

modelsub = pe.Node(util.Merge(len(model_output.outputs.__dict__)),
                   name="modelsinksub")

flutil.get_output_substitutions(model, model_output, modelsub)

# Model node substitutions
# NOTE: It would be nice if this were more intuitive, but I haven't
# figured out a good way.  Have to hardcode the node names for now.
model_mapnodes = ["featmodel", "sliceresidual", "modelestimate", "slicestats"]

modelsinksubs = flutil.get_mapnode_substitutions(exp.nruns, model_mapnodes)

modelsinksubs.append(("_contrast_","stats/"))

flutil.set_substitutions(model, modelsink, modelsub, modelsinksubs)

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
    ])

# Connect inputs
flutil.connect_inputs(model, modelsource, model_input)

# Set up the subject containers
flutil.subject_container(model, subjectsource, modelsink)

# Connect the heuristic outputs to the datainks
flutil.sink_outputs(model, model_output, modelsink, "model.%s"%spacedict[args.space])

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


# Registration
# ============

# Get registration input and output
registration, reg_input, reg_output = get_registration_workflow(surface=False)

imagesource = pe.Node(util.IdentityInterface(fields=["image"]),
                      iterables=("image", ["cope", "varcope"]),
                      name="imagesource")

# Registration datasource node
regsource = pe.Node(nio.DataGrabber(infields=["subject_id",
                                              "contrast",
                                              "image"],
                                    outfields=["warpfield",
                                               "fsl_affine",
                                               "source_images"],
                                    base_directory=analysis_dir,
                                    sort_filelist=True),
                    name="regsource")

regsource.inputs.template = "%s/%s/run_?/%s"
regsource.inputs.field_template = dict(
    warpfield=pjoin(data_dir, "%s/normalization/warpfield.nii.gz"),
    source_images="%s/model/%s/run_?/%s%d.nii.gz")
regsource.inputs.template_args = dict(fsl_affine=[["subject_id", "preproc", "func2anat_flirt.mat"]],
                                      source_images=[["subject_id", spacedict[args.space], "image", "contrast"]],
                                      warpfield=[["subject_id"]])


# Registration Datasink nodes
regsink = pe.Node(nio.DataSink(base_directory=analysis_dir),
                  name="regsink")



# Registration connections
registration.connect([
    (subjectsource, regsource,     [("subject_id", "subject_id")]),
    (imagesource,   regsource,     [("image", "image")]),
    (connames,      regsource,     [(("contrast", lambda x: [i+1 for i in get_contrast_idx(x)]), "contrast")]),
    (regsource,     reg_input,     [("warpfield", "warpfield"),
                                    ("source_images", "vol_source"),
                                    ("fsl_affine", "fsl_affine")]),
    ])

# Registration filename substitutions
regsinksub = pe.Node(util.Merge(len(reg_output.outputs.__dict__)),
                     name="regsinksub")

flutil.get_output_substitutions(registration, reg_output, regsinksub)

# Registration node substitutions
reg_mapnodes = ["applywarp", "surfproject", "surftransform",
                "convertnormsurf","smoothnativesurf", "smoothnormsurf"]

regsinksubs = flutil.get_mapnode_substitutions(exp.nruns, reg_mapnodes)
regsinksubs.extend([("_image_cope", ""),
                    ("_image_varcope", ""),
                    ("_contrast_", "")])

def set_reg_subs(contrast):
    
    for c in get_contrast_idx(contrast):
        c += 1
        regsinksubs.extend([("cope%d"%c, "cope"), ("varcope%d"%c, "varcope")])
    return regsinksubs

registration.connect(connames, ("contrast", set_reg_subs), regsink, "substitutions")

# Set up the subject containers
flutil.subject_container(registration, subjectsource, regsink)

# Connect the heuristic outputs to the datainks
flutil.sink_outputs(registration, reg_output, regsink, "registration")

# Registration working output
registration.base_dir = working_dir

# Set crashdumps
flutil.archive_crashdumps(registration)


# Across-run Fixed Effects
# ========================

# Get the workflow and nodes
fixed_fx, ffx_input, ffx_output = get_fixedfx_workflow()

# Fixed fx datasource
ffxsource = pe.Node(nio.DataGrabber(infields=["subject_id",
                                              "contrast"],
                                    outfields=["cope",
                                               "varcope",
                                               "dof_file"],
                                    base_directory=analysis_dir,
                                    sort_filelist=True),
                    name="ffxsource")

ffxsource.inputs.template = "%s/registration/%s/run_?/%s_warp.nii.gz"
ffxsource.inputs.field_template = dict(dof_file="%s/model/smoothed/run_?/dof")
ffxsource.inputs.template_args = dict(cope=[["subject_id", "contrast", "cope"]],
                                      varcope=[["subject_id", "contrast", "varcope"]],
                                      dof_file=[["subject_id"]])

# Connect inputs
fixed_fx.connect([
    (subjectsource, ffxsource, [("subject_id", "subject_id")]),
    (connames,      ffxsource, [("contrast", "contrast")]),
    ])
flutil.connect_inputs(fixed_fx, ffxsource, ffx_input, makelist = ["cope", "varcope", "dof_file"])

# Fixed effects Datasink nodes
ffxsink = pe.Node(nio.DataSink(base_directory=analysis_dir),
                  name="ffxsink")

# Fixed effects filename substitutions
ffxsinksub = pe.Node(util.Merge(len(ffx_output.outputs.__dict__)),
                     name="ffxsinksub")

flutil.get_output_substitutions(fixed_fx, ffx_output, ffxsinksub)

# Fixed effects node substitutions
flutil.set_substitutions(fixed_fx, ffxsink, ffxsinksub, [("_contrast_", ""), ("_hemi_","")])

# Set up the subject containers
flutil.subject_container(fixed_fx, subjectsource, ffxsink)

# Connect the heuristic outputs to the datainks
flutil.sink_outputs(fixed_fx, ffx_output, ffxsink, "fixed_fx.%s"%args.space)

# Fixed effects working output
fixed_fx.base_dir = working_dir

# Set crashdumps
flutil.archive_crashdumps(fixed_fx)

if __file__.endswith("fluid_fmri.py"):
    report_script = "python /mindhive/gablab/fluid/Nipype_Code/reporting/build_report.py "
else:
    report_script = "python /dev/null "

def report():
    print "Building report"
    os.system(report_script + " ".join(subject_list))
        
def workflow_runner(wf, stem):
    if any([arg for arg in args.workflows if arg.startswith(stem)]) or "all" in args.workflows:
        wf.run(inseries=args.inseries)
        report()

if __name__ == "__main__":
    if args.run:
        workflow_runner(preproc, "preproc")
        workflow_runner(model, "model") 
        workflow_runner(registration, "reg")
        workflow_runner(fixed_fx, "ffx")        
