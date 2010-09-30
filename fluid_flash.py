"""
NiPype module for the analysis of multispectral FLASH images
"""

import os
import sys
import argparse
from glob import glob
from datetime import datetime
import numpy as np

import nipype.pipeline.engine as pe

import nipype.interfaces.io as nio
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.interfaces.freesurfer as fs

parser = argparse.ArgumentParser(description="FLASH processing stream for GFluid project.")


working_base = "/mindhive/gablab/fluid/Analysis/NiPype/workingdir/flash"
data_dir = "/mindhive/gablab/fluid/Data"
analysis_dir = "/mindhive/gablab/fluid/Analysis/NiPype/flash"

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

# Subject Info Source
# -------------------

sidsource = pe.Node(util.IdentityInterface(fields=["sid"]),name="subjsource")
sidsource.iterables = ("sid", ["gf20"])

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

convert_pipe.connect(sidsource, "sid", dcmgrabber, "sid")
convert_pipe.connect(
    sidsource, ("sid", lambda x: [t[1] for t in find_flash(x)]), dcmgrabber, "dcm")
convert_pipe.connect(
    sidsource, ("sid", lambda x: [t[4] for t in find_flash(x)]), dcmgrabber, "type")

# Convert each FLASH echo to mgz format
convertdcm = pe.MapNode(fs.MRIConvert(out_type="mgz"),
                        name="convertdcm",iterfield=["args","in_file"])

convert_pipe.connect(dcmgrabber, "dicom", convertdcm, "in_file")
convert_pipe.connect(
    sidsource, ("sid", lambda x: ["-nth %d"%t[3] for t in find_flash(x)]), convertdcm, "args")

# Save the converted FLASH images in the data directory
sinkflash = pe.Node(nio.DataSink(base_directory=data_dir),name="sinkflash")

convert_pipe.connect([(sidsource, sinkflash, [("sid", "container"),
                                             (("sid", lambda x: "_sid_"+x), "strip_dir"),
                                             (("sid", get_dcm_substitutes), "substitutions")]),
                     (convertdcm, sinkflash, [("out_file", "flash.@mgz_files")])
                     ])                              

# Registration Pipeline
# ---------------------

reg_pipe = pe.Workflow(name="registration_pipe",
                       base_dir=working_base)

# Provide the set of angles for our FLASH acquisistions
anglesource = pe.Node(util.IdentityInterface(fields=["alpha"]),name="anglesource")
anglesource.iterables = ("alpha", [5, 20, 30])

# Collect FLASH images in compressed mgh format from the Data directory
flashgrabber = pe.Node(nio.DataGrabber(infields=["alpha","sid"],outfields=["flash_files"]),
                          name="flashgrabber")

flashgrabber.inputs.base_directory = os.path.join(data_dir)
flashgrabber.inputs.template = "%s/flash/flash_%02d-?.mgz"
flashgrabber.inputs.template_args = dict(flash_files=[["sid", "alpha"]])

reg_pipe.connect([(sidsource, flashgrabber, [("sid", "sid")]),
                  (anglesource, flashgrabber, [("alpha", "alpha")])
                  ])

# Get a mean image of the 8 FLASH echos
meanflash = pe.Node(fs.Concatenate(stats="mean"), name="meanflash")

reg_pipe.connect(flashgrabber, "flash_files", meanflash, "in_files")

# Skullstrip the mean images
stripflash = pe.Node(fsl.BET(), name="stripflash")

reg_pipe.connect(meanflash, "concatenated_file", stripflash, "in_file")

# Register each mean image to the Freesurfer structural
coregister = pe.Node(fs.BBRegister(init="fsl"),
                     name="coregister")

# Create a dictionary to associate expected image contrast with flip angle
angle2contrast = {5:"t2",20:"t1",30:"t1"}

reg_pipe.connect([
    (sidsource, coregister, [("sid", "subject_id")]),
    (anglesource, coregister, [(("alpha", lambda x: angle2contrast[x]), "contrast_type")]),
    (stripflash, coregister, [("out_file", "source_file")])
                  ])

# Apply the transformation to each echo image
applyxfm = pe.MapNode(fs.ApplyVolTransform(fs_target=True),
                      iterfield=["source_file"],
                      name="applyxfm")

reg_pipe.connect([(coregister, applyxfm, [("out_reg_file", "reg_file")]),
                  (flashgrabber, applyxfm, [("flash_files", "source_file")])
                  ])

# Sink the transformed images
xfmsink = pe.Node(nio.DataSink(base_directory=analysis_dir),name="xfmsink")
xfmsink.inputs.parameterization = False
xfmsink.inputs.substitutions = [("_warped","")]
reg_pipe.connect([(sidsource, xfmsink, [("sid", "container"),
                                        (("sid", lambda x: "_sid_" + x), "strip_dir")]),
                  (applyxfm, xfmsink, [("transformed_file", "coregistration.@flash_files")])
                  ])


# Parameter Estimation Pipeline
# -----------------------------

fit_pipe = pe.Workflow(name="parameter_estimation_pipe",base_dir=working_base)

# Grab all of the coregistered FLASH images from the analysis directory
regflashgrabber = pe.Node(nio.DataGrabber(infields=["sid"],outfields=["flash_files"]),
                          name="regflashgrabber")
regflashgrabber.inputs.base_directory = analysis_dir
regflashgrabber.inputs.template = "%s/coregistration/flash_??-?.mgz"
regflashgrabber.inputs.template_args = dict(flash_files=[["sid"]])

fit_pipe.connect(sidsource, "sid", regflashgrabber, "sid")

# Estimate tissue parameters
fitparams = pe.Node(fs.FitMSParams(), name="fitparams")

fit_pipe.connect(regflashgrabber, "flash_files", fitparams, "in_files")

# Convert the T1 volume to nifti
convert2nii = pe.Node(fs.MRIConvert(out_type="niigz"), name="convert2nii")

fit_pipe.connect(fitparams, "t1_image", convert2nii, "in_file")

# Get the R1 volume by taking the reciprocal of the T1 volume
invert = pe.Node(fsl.ImageMaths(op_string="-recip",suffix="recip"), name="invert_t1")

fit_pipe.connect(convert2nii, "out_file", invert, "in_file")

# Sink the T1 and R1 volumes in the analysis directory
paramsink = pe.Node(nio.DataSink(base_directory=analysis_dir),name="paramsink")
paramsink.inputs.parameterization = False
paramsink.inputs.substitutions = [("T1_out_recip","R1")]
fit_pipe.connect([(sidsource, paramsink, [("sid", "container"),
                                          (("sid", lambda x: "_sid_"+x), "strip_dir")]),
                  (fitparams, paramsink, [("t1_image", "tissue_parameters.@T1_vol")]),
                  (invert, paramsinl, [("out_file", "tissue_parameters.@R1_vol")])
                  ])

# Set the crashdump directory                           
datestamp = str(datetime.now())[:10]
codepath = os.path.split(os.path.abspath(__file__))[0]
crashdir = os.path.abspath("%s/crashdumps/%s" % (codepath, datestamp))
if not os.path.isdir(crashdir):    
    os.makedirs(crashdir)
convert_pipe.config = dict(crashdump_dir=crashdir) 
reg_pipe.config = dict(crashdump_dir=crashdir)
fit_pipe.config = dict(crashdump_dir=crashdir)

