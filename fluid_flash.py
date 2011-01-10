#! /usr/bin/env python
"""
Nipype module for the analysis of multispectral FLASH images

"""
import os
import argparse
from glob import glob
import numpy as np

import nipype.pipeline.engine as pe

import nipype.interfaces.io as nio
import nipype.interfaces.fsl as fsl
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
parser.add_argument("-inseries",action="store_true",
                    help="force running in series")
args = parser.parse_args()

if args.workflows == ["all"]:
    args.workflows = ["cvt", "reg", "fit"]

if args.subjects is None:
    subject_list = []
else:
    subject_list = args.subjects

# Hard code some paths
# --------------------
working_base = "/mindhive/gablab/fluid/Analysis/Nipype/workingdir/flash"
data_dir = "/mindhive/gablab/fluid/Data"
analysis_dir = "/mindhive/gablab/fluid/Analysis/Nipype/flash"

# Functions to control dicom file sourcing
# ----------------------------------------

def get_info_file(sid):
    """Get a dicom info file from one of the two places it could be."""
    try:
        path = glob(os.path.join(data_dir,"%s/unpack/struct-dicominfo.txt"%sid))[0]
        type = "struct"
        return path, type
    except IndexError:
        path = glob(os.path.join(data_dir,"%s/unpack/full-dicominfo.txt"%sid))[0]
        type = "full"
        return path, type

def find_flash(sid):
    """Parse a dicom info file to find flash run ids."""
    infofile, type = get_info_file(sid)
    seqinfo = np.genfromtxt(infofile,dtype=object)
    info = []
    for l in seqinfo:
        seqn = l[2]
        name=l[12]
        dcmfile=l[1]
        t = l[9]
        if name.startswith("gre_mgh_multiecho") and (t=="8"):
            alpha = name[18:-12]
            for echo in range(8):
                info.append((seqn,dcmfile,alpha,echo,type))
    return info

def get_dcm_substitutes(sid):
    """Rename the output of a the dicom conversion pipeline sensibly."""
    infolist = find_flash(sid)
    subs = []
    for i, info in enumerate(infolist):
        fname = info[1][:-4]
        alpha = int(info[2])
        echo = int(info[3])
        subs.append(("_convertdcm%d/%s_out.mgz"%(i,fname),"flash_%02d-%d.mgz"%(alpha,echo)))
    return subs

def get_angle_groups(subject_list):
    """Figure out the flip angle groups for each subject."""
    angles = {(5,20,30):[],(5,30):[]}
    for sid in subject_list:
        flashfiles = glob(os.path.join(data_dir, sid, "flash", "flash_??-?.mgz"))
        if len(flashfiles) == 24:
            angles[(5,20,30)].append(sid)
        else:
            angles[(5,30)].append(sid)
    return angles

# Subject Info Source
# -------------------

sidsource = pe.Node(util.IdentityInterface(fields=["sid"]),name="subjsource")
sidsource.iterables = ("sid", subject_list)

# Conversion Pipeline
# -------------------

convert_pipe = pe.Workflow(name="convert_pipe",
                           base_dir=working_base)

# Get the dicom files from the data directory
dcmgrabber = pe.MapNode(nio.DataGrabber(infields=["sid","type","dcm"],outfields=["dicom"]),
                        iterfield=["type","dcm"],
                        name="dcmgrabber")
dcmgrabber.inputs.base_directory="/mindhive/gablab/fluid/Data"
dcmgrabber.inputs.template="%s/dicom/%s/%s"
dcmgrabber.inputs.template_args=dict(dicom=[["sid","type","dcm"]])


# Convert each FLASH echo to mgz format
convertdcm = pe.MapNode(fs.MRIConvert(out_type="mgz"),
                        name="convertdcm",iterfield=["args","in_file"])

# Save the converted FLASH images in the data directory
sinkflash = pe.Node(nio.DataSink(base_directory=data_dir),name="sinkflash")

convert_pipe.connect([
    (sidsource,   dcmgrabber,  [("sid", "sid"),
                                (("sid", lambda x: [t[1] for t in find_flash(x)]), "dcm"),
                                (("sid", lambda x: [t[4] for t in find_flash(x)]), "type")]),
    (dcmgrabber,  convertdcm,  [("dicom", "in_file")]),
    (sidsource,   convertdcm,  [(("sid", lambda x:["-nth %d"%t[3] for t in find_flash(x)]), "args")]),
    (sidsource,   sinkflash,   [("sid", "container"),
                                (("sid", lambda x: "_sid_"+x), "strip_dir"),
                                (("sid", get_dcm_substitutes), "substitutions")]),
    (convertdcm,  sinkflash,   [("out_file", "flash.@mgz_files")]),
    ])                              

# Registration Pipeline
# ---------------------

reg_pipe = pe.Workflow(name="registration_pipe",
                       base_dir=working_base)

# Provide the set of angles for our FLASH acquisistions
anglesource = pe.Node(util.IdentityInterface(fields=["alpha"]),name="anglesource")

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

# Convert the RMS image to nifti so BET can read it
cvt2nii = pe.Node(fs.MRIConvert(out_type="niigz"), name="cvt2nii")

# Skullstrip the RMS images
stripflash = pe.Node(fsl.BET(), name="stripflash")

# Register each RMS image to the Freesurfer structural
coregister = pe.Node(fs.BBRegister(init="fsl"),
                     name="coregister")

# Create a dictionary to associate expected image contrast with flip angle
angle2contrast = {5:"t2",20:"t1",30:"t1"}

# Apply the transformation to each echo image
applyxfm = pe.MapNode(fs.ApplyVolTransform(fs_target=True),
                      iterfield=["source_file"],
                      name="applyxfm")

# Sink the transformed images
xfmsink = pe.Node(nio.DataSink(base_directory=analysis_dir),name="xfmsink")
xfmsink.inputs.parameterization = False
xfmsink.inputs.substitutions = [("_warped","")]


reg_pipe.connect([
    (sidsource,    flashgrabber, [("sid", "sid")]),
    (anglesource,  flashgrabber, [("alpha", "alpha")]),
    (sidsource,    rmsgrabber,   [("sid", "sid")]),
    (anglesource,  rmsgrabber,   [("alpha", "alpha")]),
    (rmsgrabber,   cvt2nii,      [("flash_rms", "in_file")]),
    (cvt2nii,      stripflash,   [("out_file", "in_file")]),
    (sidsource,    coregister,   [("sid", "subject_id")]),
    (anglesource,  coregister,   [(("alpha", lambda x: angle2contrast[x]), "contrast_type")]),
    (stripflash,   coregister,   [("out_file", "source_file")]),
    (coregister,   applyxfm,     [("out_reg_file", "reg_file")]),
    (flashgrabber, applyxfm,     [("flash_echos", "source_file")]),
    (sidsource,    xfmsink,      [("sid", "container"),
                                  (("sid", lambda x: "_sid_" + x), "strip_dir")]),
    (applyxfm,     xfmsink,      [("transformed_file", "coregistration.@flash_files")]),
    ])


# Parameter Estimation Pipeline
# -----------------------------

fit_pipe = pe.Workflow(name="parameter_estimation_pipe",base_dir=working_base)

# Grab all of the coregistered FLASH images from the analysis directory
regflashgrabber = pe.Node(nio.DataGrabber(infields=["sid"],outfields=["flash_echos"]),
                          name="regflashgrabber")
regflashgrabber.inputs.base_directory = analysis_dir
regflashgrabber.inputs.template = "%s/coregistration/flash_??-?.mgz"
regflashgrabber.inputs.template_args = dict(flash_echos=[["sid"]])

# Estimate tissue parameters
fitparams = pe.Node(fs.FitMSParams(), name="fitparams")

hemisource = pe.Node(util.IdentityInterface(fields=["hemi"]),name="hemisource")
hemisource.iterables = ("hemi", ["lh","rh"])

# Sample the T1 volume onto the cortical surface
t1surf = pe.Node(fs.SampleToSurface(reg_header=True,
                                    sampling_range=.5,
                                    sampling_units="frac",
                                    cortex_mask=True),
                 name="t1surf")
t1surf.inputs.sampling_method="point"

#Take screenshots of the T1 parameters on the surface
surfshots = pe.Node(fs.SurfaceScreenshots(surface="inflated",
                                          show_color_scale=True,
                                          show_gray_curv=True,
                                          six_images=True,
                                          overlay_range=(950,2500)),
                    name="surfshots")

# Sink the T1 volumes and surfaces into the analysis directory
paramsink = pe.Node(nio.DataSink(base_directory=analysis_dir),name="paramsink")
paramsink.inputs.parameterization = False
fit_pipe.connect([
    (sidsource,        regflashgrabber, [("sid", "sid")]), 
    (regflashgrabber,  fitparams,       [("flash_echos", "in_files")]),
    (fitparams,        t1surf,          [("t1_image", "source_file")]),
    (sidsource,        t1surf,          [("sid", "subject_id")]),
    (hemisource,       t1surf,          [("hemi", "hemi")]),
    (sidsource,        surfshots,       [("sid", "subject")]),
    (hemisource,       surfshots,       [("hemi", "hemi")]),
    (t1surf,           surfshots,       [("out_file", "overlay")]),
    (sidsource,        paramsink,       [("sid", "container"),
                                         (("sid", lambda x: "_sid_"+x), "strip_dir")]),
    (fitparams,        paramsink,       [("t1_image", "tissue_parameters.@T1_vol")]),
    (t1surf,           paramsink,       [("out_file", "tissue_parameters.@T1_surf")]),
    (surfshots,        paramsink,       [("screenshots", "screenshots.@images")]),
    ])



# Set the crashdump directory                           
for pipe in [convert_pipe, reg_pipe, fit_pipe]:
    flutil.archive_crashdumps(pipe)

def run_reg_pipe():
    angles = get_angle_groups(subject_list)
    for anglegroup in angles:
        anglesource.iterables = ("alpha", anglegroup)
        sidsource.iterables = ("sid", angles[anglegroup])
        reg_pipe.run()
    sidsource.iterables = ("sid", subject_list)

if __name__ == "__main__":
    if "cvt" in args.workflows:
        convert_pipe.run()
    if "reg" in args.workflows:
        run_reg_pipe()
    if "fit" in args.workflows:
        fit_pipe.run()
