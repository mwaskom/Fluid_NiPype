#! /usr/bin/env python
"""
Main interface for GFluid fMRI Nipype code.
"""
import os
import sys
import argparse
import inspect
from warnings import warn
from copy import deepcopy

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
from nipype.interfaces.base import Bunch

from workflows.preproc import preproc
from workflows.registration import registration
from workflows.fsl_model import volume_model
from workflows.surface_projection import surfproj
from workflows.fixed_fx_vol import fixed_fx as fixed_fx_volume
from workflows.fixed_fx_surf import fixed_fx as fixed_fx_surface

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
# Look for paradgim_experiment.py in an experiments package
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
preprocsinksub = pe.Node(util.Merge(len(preproc_output.outputs.__dict__)),
                         name="preprocsinksub")

flutil.get_output_substitutions(preproc, preproc_output, preprocsinksub)

preprocreportsub = pe.Node(util.Merge(len(preproc_report.outputs.__dict__)),
                              name="preprocreportsub")

flutil.get_output_substitutions(preproc, preproc_report, preprocreportsub)

# Preproc node substitutions
# NOTE: It would be nice if this were more intuitive, but I haven't
# figured out a good way.  Have to hardcode the node names for now.
preproc_mapnodes = ["art", "dilatemask", "highpass", "masksmoothfunc", 
                    "extractref", "meanfunc2", "realign"]
preprocsinknodesubs = flutil.get_mapnode_substitutions(exp.nruns, preproc_mapnodes)

preproc_report_mapnodes = []
for plot in ["displacement", "rotation", "translation"]:
    preproc_report_mapnodes.append("plot%s"%plot)
for img in ["example", "mean"]:
    preproc_report_mapnodes.append("%sslice"%img)
preprocreportnodesubs = flutil.get_mapnode_substitutions(exp.nruns, preproc_report_mapnodes)


flutil.set_substitutions(preproc, preprocsink, preprocsinksub, preprocsinknodesubs)

flutil.set_substitutions(preproc, preprocreport, preprocreportsub, preprocreportnodesubs)

# Preproc connections
preproc.connect(subjectsource, "subject_id", preprocsource, "subject_id")

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
                                               "functional_mask", "example_func",
                                               "vol_timeseries", "surf_timeseries"],
                                    base_directory = project_dir,
                                    sort_filelist=True),
                    name="regsource")

regsource.inputs.template = "Analysis/NiPype/" + args.paradigm + "/%s/preproc/run_?/%s.nii.gz"
regsource.inputs.field_template = dict(warpfield="Data/%s/normalization/%s.nii.gz")
regsource.inputs.template_args = dict(mean_func=[["subject_id", "mean_func"]],
                                      example_func=[["subject_id", "example_func"]],
                                      functional_mask=[["subject_id", "functional_mask"]],
                                      vol_timeseries=[["subject_id", "smoothed_timeseries"]],
                                      surf_timeseries=[["subject_id", "unsmoothed_timeseries"]],
                                      warpfield=[["subject_id","warpfield"]])

# Registration Datasink nodes
regsink = pe.Node(nio.DataSink(base_directory=analysis_dir),
                  name="regsink")

regreport = pe.Node(nio.DataSink(base_directory=report_dir),
                    name="regreport")

# Registration connections
registration.connect([
    (subjectsource, regsource, [("subject_id", "subject_id")]),
    (subjectsource, reg_input, [("subject_id", "subject_id")]),
    ])

# Assign a smoothing value for the registration workflow
reg_input.inputs.smooth_fwhm = 0.
warn("\nWARNING:\n\nYou haven't decided on a surface smoothing kernel yet")

# Registration filename substitutions
regsinksub = pe.Node(util.Merge(len(reg_output.outputs.__dict__)),
                     name="regsinksub")

flutil.get_output_substitutions(registration, reg_output, regsinksub)

regreportsub = pe.Node(util.Merge(len(reg_report.outputs.__dict__)),
                       name="regreportsub")

flutil.get_output_substitutions(registration, reg_report, regreportsub)

# Registration node substitutions
reg_mapnodes = ["func2anat","warptimeseries","warpexample","warpmask","meanwarp",
                "convertnormsurf","smoothnativesurf"]
regsinknodesubs = flutil.get_mapnode_substitutions(exp.nruns, reg_mapnodes)

flutil.set_substitutions(registration, regsink, regsinksub, regsinknodesubs)

reg_report_mapnodes = ["exfuncwarppng", "func2anat", "func2anatpng"]
regreportnodesubs = flutil.get_mapnode_substitutions(exp.nruns, reg_report_mapnodes)

flutil.set_substitutions(registration, regreport, regreportsub, regreportnodesubs)

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
                            amplitudes=None,
                            tmod=None,
                            pmod=None,
                            regressor_names=None,
                            regressors=None))
    return output


# Set up the model workflow for each space 
model = volume_model.clone("%s_model"%args.space)

# Get model input and output
model_input = model.get_node("inputspec")
model_output = model.get_node("outputspec")
model_report = model.get_node("report")

# Set up the date sources for each space
modelsource = pe.Node(nio.DataGrabber(infields=["subject_id"],
                                      outfields=["timeseries",
                                                 "outlier_files",
                                                 "overlay_background",
                                                 "realignment_parameters"],
                                      base_directory = project_dir,
                                      sort_filelist=True),
                      name="modelsource")

modelsource.inputs.template = os.path.join(
    "Analysis/NiPype",args.paradigm ,"%s/%s/run_?/%s")

modelsource.inputs.template_args = dict(
    outlier_files=[["subject_id", "preproc", "outlier_files.txt"]],
    realignment_parameters=[["subject_id", "preproc", "realignment_parameters.par"]])

# Set the template and args for the timeseries in each space
if args.space == "volume":
    modelsource.inputs.template_args["timeseries"] = [["subject_id", "registration", "warped_timeseries.nii.gz"]]
    modelsource.inputs.field_template = dict(overlay_background=fsl.Info().standard_image("avg152T1.nii.gz"))
    modelsource.inputs.template_args["overlay_background"] = [[]]
elif args.space == "surface":
    modelsource.inputs.template_args["timeseries"] = [["subject_id", "preproc", "unsmoothed_timeseries.nii.gz"]]
    modelsource.inputs.template_args["overlay_background"] = [["subject_id", "preproc", "example_func.nii.gz"]]
else:
    sys.exit("Model space %s not recognized"%args.space)



# Model Datasink nodes
modelsink = pe.Node(nio.DataSink(base_directory=analysis_dir),
                  name="modelsink")

modelreport = pe.Node(nio.DataSink(base_directory=report_dir),
                    name="modelreport")

modelreportsub = pe.Node(util.Merge(len(model_report.outputs.__dict__)),
                       name="modelreportsub")

flutil.get_output_substitutions(model, model_report, modelreportsub)

# Model node substitutions
# NOTE: It would be nice if this were more intuitive, but I haven't
# figured out a good way.  Have to hardcode the node names for now.
modelsinknodesubs = []
for r in range(exp.nruns):
    for node in ["modelestimate"]:
        modelsinknodesubs.append(("_%s%d"%(node, r), "run_%d"%(r+1)))
modelsink.inputs.substitutions = modelsinknodesubs

modelreport_mapnodes = ["slicestats", "modelgen", "sliceresidual"]

modelreportnodesubs = flutil.get_mapnode_substitutions(exp.nruns, modelreport_mapnodes)

modelreportnodesubs.append(("_contrast_","stats/"))


flutil.set_substitutions(model, modelreport, modelreportsub, modelreportnodesubs)

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
flutil.connect_inputs(model, modelsource, model_input, makelist=["overlay_background"], listlength=exp.nruns)

# Set up the subject containers
flutil.subject_container(model, subjectsource, modelsink)
flutil.subject_container(model, subjectsource, modelreport)

# Connect the heuristic outputs to the datainks
flutil.sink_outputs(model, model_output, modelsink, "model.%s"%args.space)
flutil.sink_outputs(model, model_report, modelreport, "model.%s"%args.space)

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

# Sampling to Surface
# -------------------

# Get sampling input and output
surfproj_input = surfproj.get_node("inputspec")
surfproj_output = surfproj.get_node("outputspec")

# Surface sampling datasource
surfsource = pe.Node(nio.DataGrabber(infields=["subject_id"],
                                     outfields=["reg_matrix",
                                                "cope",
                                                "varcope"],
                                     base_directory=project_dir,
                                     sort_filelist=True),
                     name="surfsource")

if args.space == "surface":
    surfsource.inputs.template = "Analysis/NiPype/" + args.paradigm + "/%s/%s/run_?/%s%d.nii.gz"
    surfsource.inputs.field_template = dict(
        reg_matrix= "Analysis/NiPype/" + args.paradigm + "/%s/registration/run_?/register.dat")
    surfsource.inputs.template_args = dict(cope = [["subject_id", "model/surface", "cope", "contrast"]],
                                           varcope = [["subject_id", "model/surface", "varcope", "contrast"]],
                                           reg_matrix = [["subject_id"]])

# Connect contrast names to datasource, which wants a 1-based index
def get_con_img_idx(contrast_name):
    """Return the index corresponding to the name of a contrast."""
    return [i+1 for i, data in enumerate(exp.contrasts) if data[0] == contrast_name][0]

surfproj.connect(connames, ("contrast", get_con_img_idx), surfsource, "contrast")

# Connect the inputs
flutil.connect_inputs(surfproj, surfsource, surfproj_input, makelist = ["volume_file"])

surfproj.connect([
    (subjectsource, surfsource,   [("subject_id", "subject_id")]),
    (subjectsource, surfproj_input, [("subject_id", "subject_id")]),
    ])

# Set the smoothing value
surfproj_input.inputs.smooth_fwhm = 5.

# Surf Projection Datasink
surfprojsink = pe.Node(nio.DataSink(base_directory=analysis_dir),
                       name="surfprojsink")

# Filename substitutions
surfprojsinksub = pe.Node(util.Merge(len(surfproj_output.outputs.__dict__)),
                          name="surfprojsinksub")

flutil.get_output_substitutions(surfproj, surfproj_output, surfprojsinksub)

flutil.sink_outputs(surfproj, surfproj_output, surfprojsink, "surfprojection")

# Node substitutions
surfproj_mapnodes = ["copesmoother","varcopesmoother"]
surfprojsinknodesubs = flutil.get_mapnode_substitutions(exp.nruns, surfproj_mapnodes)
surfprojsinknodesubs.extend([("_contrast_",""),("_hemi_","")])

flutil.set_substitutions(surfproj, surfprojsink, surfprojsinksub, surfprojsinknodesubs)

# Set up the subject containers
flutil.subject_container(surfproj, subjectsource, surfprojsink)

# Working output
surfproj.base_dir = working_dir

# Set crashdumps
flutil.archive_crashdumps(surfproj)


# Across-run Fixed Effects
# ------------------------

# Figure out which fixed_fx workflow we want
fixed_fx = locals()["fixed_fx_%s"%args.space]

# "Use" the fixed_fx workflows so pyflake shuts up
_ = fixed_fx_volume
_ = fixed_fx_surface

# Get fixed_fx input and output
ffx_input = fixed_fx.get_node("inputspec")
ffx_output = fixed_fx.get_node("outputspec")
ffx_report = fixed_fx.get_node("report")

# Fixed effects datasource node
if args.space == "volume":
    ffxsource = pe.Node(nio.DataGrabber(infields=["subject_id","contrast"],
                                        outfields=["cope","varcope"],
                                        base_directory = project_dir,
                                        sort_filelist=True),
                        name="ffxsource")

    ffxsource.inputs.template = "Analysis/NiPype/" + args.paradigm + "/%s/model/volume/run_?/%s%d.nii.gz"
    ffxsource.inputs.template_args = dict(cope = [["subject_id", "cope", "contrast"]],
                                          varcope = [["subject_id", "varcope", "contrast"]])
elif args.space == "surface":
    ffxsource = pe.Node(nio.DataGrabber(infields=["subject_id","contrast","hemi"],
                                        outfields=["cope"],
                                        base_directory = project_dir,
                                        sort_filelist=True),
                        name="ffxsource")

    ffxsource.inputs.template = "Analysis/NiPype/" + args.paradigm + "/%s/surfprojection/%s/%s/run_?/%s.nii.gz"
    ffxsource.inputs.template_args = dict(cope = [["subject_id", "contrast", "hemi", "cope"]],
                                          varcope = [["subject_id", "contrast", "hemi", "varcope"]])

    hemisource = pe.Node(util.IdentityInterface(fields=["hemi"]),
                         iterables=("hemi",["lh","rh"]),
                         name="hemisource")
    fixed_fx.connect(hemisource, "hemi", ffxsource, "hemi")

# Fixed effects Datasink nodes
ffxsink = pe.Node(nio.DataSink(base_directory=analysis_dir),
                  name="ffxsink")

ffxreport = pe.Node(nio.DataSink(base_directory=report_dir),
                    name="ffxreport")

# Fixed effects connections
fixed_fx.connect(subjectsource, "subject_id", ffxsource, "subject_id")
if args.space == "surface":
    fixed_fx.connect(hemisource, "hemi", ffx_input, "hemi")

# Connect contrast names to datasource
fixed_fx.connect(connames, "contrast", ffxsource, "contrast")

# Fixed effects filename substitutions
ffxsinksub = pe.Node(util.Merge(len(ffx_output.outputs.__dict__)),
                     name="ffxsinksub")

flutil.get_output_substitutions(fixed_fx, ffx_output, ffxsinksub)

ffxreportsub = pe.Node(util.Merge(len(ffx_report.outputs.__dict__)),
                       name="ffxreportsub")

flutil.get_output_substitutions(fixed_fx, ffx_report, ffxreportsub)

# Fixed effects node substitutions
flutil.set_substitutions(fixed_fx, ffxsink, ffxsinksub, [("_contrast_", ""), ("_hemi_","")])

flutil.set_substitutions(fixed_fx, ffxreport, ffxreportsub, [("_contrast_", ""), ("_hemi_","")])

# Connect the inputs
flutil.connect_inputs(fixed_fx, ffxsource, ffx_input, makelist = ["cope", "varcope"])

# Set up the subject containers
flutil.subject_container(fixed_fx, subjectsource, ffxsink)
flutil.subject_container(fixed_fx, subjectsource, ffxreport)

# Connect the heuristic outputs to the datainks
flutil.sink_outputs(fixed_fx, ffx_output, ffxsink, "fixed_fx")
flutil.sink_outputs(fixed_fx, ffx_report, ffxreport, "fixed_fx")

# Fixed effects working output
fixed_fx.base_dir = working_dir

# Set crashdumps
flutil.archive_crashdumps(fixed_fx)

def report():
    if __file__ == "fluid_fmri.py":
        report_script = "python /mindhive/gablab/fluid/Nipype_Code/reporting/build_report.py "
    else:
        report_script = "python /dev/null "
    os.system(report_script + " ".join(subject_list))
        

if __name__ == "__main__":
    if args.run:
        if "preproc" in args.workflows:
            preproc.run(inseries=args.inseries)
            report()
        
        if "reg" in args.workflows:
            registration.run(inseries=args.inseries)
            report()
        
        if "model" in args.workflows:
            model.run(inseries=args.inseries)
            report()
        
        if "surfproj" in args.workflows:
            surfproj.run(inseries=args.inseries)
            report()

        if "ffx" in args.workflows:
            fixed_fx.run(inseries=args.inseries)
            report()
