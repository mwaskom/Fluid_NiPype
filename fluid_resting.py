#! /usr/bin/env python
"""
Preprocessing module for Fluid Intelligence fMRI paradigms.

Input spec node takes three inputs:
    - Timeseries (image files)
    - Highpass filter cutoff (in TRs)
    - FWHM of smoothing kernel for SUSAN (in mms)

Output spec node has two outputs:
    - Smoothed timeseries (fully preprocessed and smoothed timeseries in native space)
    - Unsmoothed timeseries (identical steps except no smoothing in the volume)
    - Example func (target volume for MCFLIRT realignment)
    - Mean func (unsmoothed mean functional)
    - Funcational mask (binary dilated brainmask in functional space)
    - Realignment parameters (text files from MCFLIRT)
    - Outlier Files (outlier text files from ART)

Reporting has two outputs:
    - Motion parameter plots (rotations and translations)
    - Motion displacement plots (absolute and relative displacement)

"""
import os
from glob import glob
import argparse
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util   
import nipype.pipeline.engine as pe         

from workflows.resting_preproc import preproc
from workflows.registration import registration

import fluid_utility as flutil

# Parse command line arguments
parser = argparse.ArgumentParser(description="Resting-state preprocessing for GFluid project")
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

project_dir = "/mindhive/gablab/fluid/"
data_dir = "/mindhive/gablab/fluid/Data"
analysis_dir = "/mindhive/gablab/fluid/Analysis/NiPype/resting"
working_dir = "/mindhive/gablab/fluid/Analysis/NiPype/workingdir/resting"
report_dir = "/mindhive/gablab/fluid/Analysis/NiPype/resting/report"

# Get all of the pre subjects as a default subject list 
if args.subjects is None:
    paths = glob(os.path.join(data_dir, "gf??"))
    paths = [p for p in paths if os.path.exists(os.path.join(p, "bold", "Resting.nii.gz"))]
    subject_list = [p.split(os.path.sep)[-1] for p in paths]
else:
    subject_list = args.subjects

print "Subjects: ", " ".join(subject_list)

if args.workflows is None:
    args.workflows = []
elif args.workflows == ["all"]:
    args.workflows = ["preproc", "reg"]

# Subject source node
# -------------------
subjectsource = pe.Node(util.IdentityInterface(fields=["subject_id"]),
                        iterables=("subject_id", subject_list),
                        name="subjectsource")

# Preprocessing
# -------------

# Get preproc input and output
preproc_input = preproc.get_node("inputspec")
preproc_output = preproc.get_node("outputspec")
preproc_report = preproc.get_node("report")

# Define a data source to grab input data
preproc_source = pe.Node(nio.DataGrabber(infields=["subject_id"],
                                        outfields=["timeseries", "anatomical", "warpfield"],
                                        base_directory=data_dir),
                        name="preprocsource")

preproc_source.inputs.template = "%s/%s/%s.nii.gz"
preproc_source.inputs.template_args = dict(timeseries = [["subject_id", "bold", "Resting"]],
                                           anatomical = [["subject_id", "normalization", "T1_fnirted"]])

                    
# Set the smoothing value
preproc_input.inputs.smooth_fwhm = 5.

# Connect the inputs
preproc.connect(subjectsource, "subject_id", preproc_source, "subject_id")
flutil.connect_inputs(preproc, preproc_source, preproc_input)

# Sink the outputs
preprocsink = pe.Node(nio.DataSink(base_directory=analysis_dir), name="datasink")
preprocreport = pe.Node(nio.DataSink(base_directory=report_dir), name="reportsink")

preprocsinksubnode = pe.Node(util.Merge(len(preproc_output.outputs.__dict__.keys())),
                             name="sinksubnode")
preprocreportsubnode = pe.Node(util.Merge(len(preproc_report.outputs.__dict__.keys())), 
                               name="reportsubnode")

# Connect the outputs to the sink
flutil.subject_container(preproc, subjectsource, preprocsink)
flutil.subject_container(preproc, subjectsource, preprocreport)
flutil.sink_outputs(preproc, preproc_output, preprocsink, "preproc")
flutil.sink_outputs(preproc, preproc_report, preprocreport, "preproc")


flutil.get_output_substitutions(preproc, preproc_output, preprocsinksubnode)
flutil.get_output_substitutions(preproc, preproc_report, preprocreportsubnode)

# Make substitutions 
preproc_mapnodes = ["art", "dilatemask", "masksmoothfunc", "maskfunc2", 
                    "extractref", "meanfunc2", "realign"]
preprocsinknodesubs = flutil.get_mapnode_substitutions(1, preproc_mapnodes)

# Change the exentsion on the realignment files so it plays nice with Conn
preprocsinknodesubs.append((".par", ".txt"))

preproc_report_mapnodes = ["realign", "art", "plotmean"]
for plot in ["displacement", "rotation", "translation"]:
    preproc_report_mapnodes.append("plot%s"%plot)
for img in ["example", "mean"]:
    preproc_report_mapnodes.append("%sslice"%img)
preprocreportnodesubs = flutil.get_mapnode_substitutions(1, preproc_report_mapnodes)


flutil.set_substitutions(preproc, preprocsink, preprocsinksubnode, preprocsinknodesubs)
flutil.set_substitutions(preproc, preprocreport, preprocreportsubnode, preprocreportnodesubs)

# Set the base directory
preproc.base_dir = os.path.join(working_dir)

# Archive crashdumps
flutil.archive_crashdumps(preproc)


# Registration
# ------------

# Get registration input and output
reg_input = registration.get_node("inputspec")
reg_output = registration.get_node("outputspec")
reg_report = registration.get_node("report")

# Get a reference to the warp node and tell it to write an uncompressed nifti
warp = registration.get_node("warptimeseries")
warp.inputs.output_type = "NIFTI"

# Registration datasource node
regsource = pe.Node(nio.DataGrabber(infields=["subject_id"],
                                    outfields=["warpfield", "mean_func",
                                               "functional_mask", "example_func",
                                               "vol_timeseries", "surf_timeseries"],
                                    base_directory = project_dir,
                                    sort_filelist=True),
                    name="regsource")

regsource.inputs.template = "Analysis/NiPype/resting/%s/preproc/run_?/%s.nii.gz"
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

# Registration filename substitutions
regsinksub = pe.Node(util.Merge(len(reg_output.outputs.__dict__.keys())),
                     name="regsinksub")

flutil.get_output_substitutions(registration, reg_output, regsinksub)

regreportsub = pe.Node(util.Merge(len(reg_report.outputs.__dict__.keys())),
                       name="regreportsub")

flutil.get_output_substitutions(registration, reg_report, regreportsub)

# Registration node substitutions
reg_mapnodes = ["func2anat","warptimeseries","warpexample","warpmask","meanwarp"]
regsinknodesubs = flutil.get_mapnode_substitutions(1, reg_mapnodes)

flutil.set_substitutions(registration, regsink, regsinksub, regsinknodesubs)

reg_report_mapnodes = ["exwarppng", "warpmaskpng", "func2anat", "func2anatpng"]
regreportnodesubs = flutil.get_mapnode_substitutions(1, reg_report_mapnodes)

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

# Run the pipeline
# ----------------

if __name__ == "__main__":
    if args.run:
        if __file__ == "fluid_resting.py":
            report_script = "python /mindhive/gablab/fluid/Nipype_Code/reporting/build_report.py "
        else:
            report_script = "python /dev/null "
        
        if "preproc" in args.workflows:
            preproc.run(inseries=args.inseries)
            os.system(report_script + " ".join(subject_list))
        
        if "reg" in args.workflows:
            registration.run(inseries=args.inseries)
            os.system(report_script + " ".join(subject_list))
