#! /usr/bin/env python
"""
Nipype module for the analysis of multispectral FLASH images

"""
import os
from os.path import join as pjoin
import argparse

import numpy as np

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
import nipype.interfaces.freesurfer as fs

import fluid_utility as flutil

# Parse the command line
# ----------------------
parser = argparse.ArgumentParser(description="FLASH processing stream for GFluid project.")
parser.add_argument("-subjects", nargs="*",
                    metavar="subject_id",
                    help="process subject(s)")
parser.add_argument("-workflows", nargs="*",
                    metavar="WF",
                    help="which workflows to run (cvt, reg, fit)")
parser.add_argument("-ipython",action="store_true")
parser.add_argument("-multiproc",action="store_true")
args = parser.parse_args()

if args.workflows == ["all"]:
    args.workflows = ["cvt", "reg", "fit"]

if args.subjects is None:
    subject_list = []
elif os.path.isfile(args.subjects[0]):
    subject_list = np.loadtxt(args.subjects[0], str).tolist()
else:
    subject_list = args.subjects

# Hard code some paths
# --------------------
project_dir =  "/mindhive/gablab/fluid/"
working_base = pjoin(project_dir, "Analysis/Nipype/workingdir/flash")
data_dir = pjoin(project_dir, "Data")
analysis_dir = pjoin(project_dir, "Analysis/Nipype/flash")

# Subject Info Source
# -------------------

sidsource = pe.Node(util.IdentityInterface(fields=["sid"]),name="subjsource")
sidsource.iterables = ("sid", subject_list)

# Registration Pipeline
# ---------------------

reg_pipe = pe.Workflow(name="registration_pipe",
                       base_dir=working_base)

# Provide the set of angles for our FLASH acquisistions
anglesource = pe.Node(util.IdentityInterface(fields=["alpha"]),name="anglesource",iterables=("alpha", [5,20,30]))

# Collect FLASH images in compressed mgh format from the Data directory
flashgrabber = pe.Node(nio.DataGrabber(infields=["alpha","sid"],
                                       outfields=["flash_echos"]),
                       name="flashgrabber")

flashgrabber.inputs.base_directory = os.path.join(data_dir)
flashgrabber.inputs.template = "%s/flash/flash_%02d-?.mgz"
flashgrabber.inputs.template_args = dict(flash_echos=[["sid", "alpha"]])


# Collect the RMS Flash images that were converted earlier
rmsgrabber = pe.Node(nio.DataGrabber(infields=["alpha","sid"],outfields=["flash_rms"]),
                     name="rmsgrabber")

rmsgrabber.inputs.base_directory = os.path.join(data_dir)
rmsgrabber.inputs.template = "%s/flash/flash_%02d-rms.mgz"
rmsgrabber.inputs.template_args = dict(flash_rms=[["sid", "alpha"]])

# Register each RMS image to the Freesurfer structural
coregister = pe.Node(fs.BBRegister(init="header"),
                     name="coregister")

# Apply the transformation to each echo image
applyxfm = pe.MapNode(fs.ApplyVolTransform(fs_target=True),
                      iterfield=["source_file"],
                      name="applyxfm")

# Get the subject's brain image
getbrain = pe.Node(nio.FreeSurferSource(subjects_dir=data_dir), name="getbrain")

# Use the brain image to skullstrip the aligned FLASH files
maskflash = pe.MapNode(fs.ApplyMask(mask_thresh=1),
                       iterfield=["in_file"],
                       name="maskflash")

# Sink the transformed images
xfmsink = pe.Node(nio.DataSink(base_directory=analysis_dir),name="xfmsink")
xfmsink.inputs.parameterization = False
xfmsink.inputs.substitutions = [("_warped_masked","")]

def angle2contrast(alpha):
    return {5:"t2",20:"t1",30:"t1"}[alpha]

reg_pipe.connect([
    (sidsource,    flashgrabber, [("sid", "sid")]),
    (anglesource,  flashgrabber, [("alpha", "alpha")]),
    (sidsource,    rmsgrabber,   [("sid", "sid")]),
    (anglesource,  rmsgrabber,   [("alpha", "alpha")]),
    (rmsgrabber,   coregister,   [("flash_rms", "source_file")]),
    (sidsource,    coregister,   [("sid", "subject_id")]),
    (anglesource,  coregister,   [(("alpha", angle2contrast), "contrast_type")]),
    (coregister,   applyxfm,     [("out_reg_file", "reg_file")]),
    (flashgrabber, applyxfm,     [("flash_echos", "source_file")]),
    (sidsource,    xfmsink,      [("sid", "container")]),
    (sidsource,    getbrain,     [("sid", "subject_id")]),
    (getbrain,     maskflash,    [("brain", "mask_file")]),
    (applyxfm,     maskflash,    [("transformed_file", "in_file")]),
    (maskflash,    xfmsink,      [("out_file", "coregistration.@flash_files")]),
    ])


# Parameter Estimation Pipeline
# -----------------------------

fit_pipe = pe.Workflow(name="parameter_estimation_pipe",base_dir=working_base)

# Grab all of the coregistered FLASH images from the analysis directory
regflashgrabber = pe.Node(nio.DataGrabber(infields=["sid"],
                                          outfields=["flash_echos"]),
                          name="regflashgrabber")
regflashgrabber.inputs.base_directory = analysis_dir
regflashgrabber.inputs.template = "%s/coregistration/flash_??-?.mgz"
regflashgrabber.inputs.template_args = dict(flash_echos=[["sid"]])

# Estimate tissue parameters
fitparams = pe.Node(fs.FitMSParams(args="-n 0 -nosynth"), name="fitparams")

hemisource = pe.Node(util.IdentityInterface(fields=["hemi"]),name="hemisource")
hemisource.iterables = ("hemi", ["lh","rh"])

# Sample the T1 volume onto the cortical surface
t1surf = pe.Node(fs.SampleToSurface(reg_header=True,
                                    sampling_range=(.1,.75,.1),
                                    sampling_units="frac",
                                    cortex_mask=True),
                 name="t1surf")
t1surf.inputs.sampling_method="average"

# Transform the native surface to fsaverage space
surfxfm = pe.Node(fs.SurfaceTransform(target_subject="fsaverage"),
                  name="surfacexfm")


# Sink the T1 volumes and surfaces into the analysis directory
paramsink = pe.Node(nio.DataSink(base_directory=analysis_dir),name="paramsink")
paramsink.inputs.parameterization = False
fit_pipe.connect([
    (sidsource,        regflashgrabber, [("sid", "sid")]), 
    (regflashgrabber,  fitparams,       [("flash_echos", "in_files")]),
    (fitparams,        t1surf,          [("t1_image", "source_file")]),
    (sidsource,        t1surf,          [("sid", "subject_id")]),
    (hemisource,       t1surf,          [("hemi", "hemi")]),
    (t1surf,           surfxfm,         [("out_file", "source_file")]),
    (sidsource,        surfxfm,         [("sid", "source_subject")]),
    (hemisource,       surfxfm,         [("hemi", "hemi")]),
    (sidsource,        paramsink,       [("sid", "container")]),
    (fitparams,        paramsink,       [("t1_image", "tissue_parameters.@T1_vol")]),
    (t1surf,           paramsink,       [("out_file", "tissue_parameters.@T1_surf")]),
    (surfxfm,          paramsink,       [("out_file", "tissue_parameters.@T1_fsaverage")]),
    ])



# Set the crashdump directory                           
for pipe in [reg_pipe, fit_pipe]:
    flutil.archive_crashdumps(pipe)

if args.ipython:
    plugin = "IPython"
elif args.multiproc:
    plugin = "MultiProc"
else:
    plugin = "Linear"


def workflow_runner(flow, stem):
    if any([a.startswith(stem) for a in args.workflows]) or args.workflows==["all"]:
        flow.run(plugin=plugin)

if __name__ == "__main__":
    workflow_runner(reg_pipe, "reg")
    workflow_runner(fit_pipe, "fit")
