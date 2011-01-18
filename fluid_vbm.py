#! /usr/bin/env python
import os, sys
from os.path import join as pjoin
import argparse

import nipype.interfaces.io as nio
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe

import fluid_utility as flutil

# Parse command line arguments
parser = argparse.ArgumentParser(description="Nipype script to run FSL VBM analysis.")

parser.add_argument("-subjects", nargs="*",
                    metavar="subject_id",
                    help="process subject(s)")
parser.add_argument("-workflows", nargs="*", metavar="<wf>",
                    help="which workflows to run")
parser.add_argument("-templatefile", metavar="<file>",
                    help="text file defining subjects to use in template")
parser.add_argument("-inseries",action="store_true",
                    help="force running in series")
args = parser.parse_args()

# Template file doesn't work yet
if args.templatefile is not None:
    sys.exit("Template file not yet implemented")

# Set up some paths
project_dir = "/mindhive/gablab/fluid"
data_dir = pjoin(project_dir, "Data")
nipype_dir = pjoin(project_dir, "Analysis/Nipype")
analysis_dir = pjoin(nipype_dir, "vbm")
working_dir = pjoin(nipype_dir, "workingdir", "vbm")

# Subject source
subjsource = pe.Node(util.IdentityInterface(fields=["subject_id"]),
                     iterables=("subject_id", args.subjects),
                     name="subjsource")

# Tissue Segmentation
# -------------------

# Grab the MPRAGEs
segsource = pe.Node(nio.DataGrabber(infields=["subject_id"],
                                    outfields=["mprage"],
                                    base_directory=data_dir,
                                    template="%s/structural/mprage.nii.gz"),
                    name="segsource")
segsource.inputs.template_args=dict(mprage=[["subject_id"]])

# Skull-strip in bias-removal mode
skullstrip = pe.Node(fsl.BET(reduce_bias=True,
                             frac=0.4),
                     name="skullstrip")

# Segment the brains
segment = pe.Node(fsl.FAST(mixel_smooth=0.3,
                           hyper=0.1),
                  name="segment")

# Sink the brains from the segmentation flow
segsink = pe.Node(nio.DataSink(base_directory=analysis_dir,
                               substitutions=[("_subject_id_",""),
                                              ("_brain", ""),
                                              ("pve_0", "csf"),
                                              ("pve_1", "gm"),
                                              ("pve_2", "wm")]),
                  name="segsink")

# Define and connect the BET workflow
seg = pe.Workflow(name="seg", base_dir=working_dir)
flutil.archive_crashdumps(seg)
seg.connect([
    (subjsource,  segsource,   [("subject_id", "subject_id")]),
    (segsource,   skullstrip,  [("mprage", "in_file")]),
    (skullstrip,  segment,     [("out_file", "in_files")]),
    (skullstrip,  segsink,     [("out_file", "@brain")]),
    (segment,     segsink,     [("partial_volume_files", "@pves")]),
    ])

# Template Creation
# -----------------

# Template subject source (could be different from
# analysis subject list, but we haven't implemented
# that yet.
tempsubjects = pe.Node(util.IdentityInterface(fields=["subject_id"]),
                       iterables = ("subject_id", args.subjects),
                       name="tempsubjects")

# Grab the grey matter partial volume images
graysource = pe.Node(nio.DataGrabber(infields=["subject_id"],
                                     outfields=["gray_image"],
                                     base_directory=analysis_dir,
                                     sort_filelist=True,
                                     template="%s/mprage_gm.nii.gz"),
                     name="graysource")
graysource.inputs.template_args=dict(gray_image=[["subject_id"]])

# Perform an initial affine transform to the gray matter tissue prior
initreg = pe.Node(fsl.FLIRT(reference=fsl.Info.standard_image("tissuepriors/avg152T1_gray.img"),
                               searchr_x=[-180,180],
                               searchr_y=[-180,180],
                               searchr_z=[-180,180]),
                     iterfield=["in_file"],
                     name="initreg")

# Sink the initial registrations
initregsink = pe.Node(nio.DataSink(base_directory=analysis_dir,
                                   substitutions=[("_subject_id_","")]),
                      name="initregsink")

# Define the intial registration flow
initregflow = pe.Workflow(name="initialreg",base_dir=working_dir)
flutil.archive_crashdumps(initregflow)
initregflow.connect([
    (tempsubjects,  graysource,  [("subject_id", "subject_id")]),
    (graysource ,   initreg,     [("gray_image", "in_file")]),
    (initreg,       initregsink, [("out_file", "@init_xfm")]),
    ])

# Grab the flirted gray matter images
initsource = pe.MapNode(nio.DataGrabber(infields=["subject_id"],
                                        outfields=["gray_image"],
                                        base_directory=analysis_dir,
                                        sort_filelist=True,
                                        template="%s/mprage_gm_flirt.nii.gz"),
                        iterfield=["subject_id"],
                        name="initsource")
initsource.inputs.template_args=dict(gray_image=[["subject_id"]])
initsource.inputs.subject_id = args.subjects

# Merge the affine-registered gray matter images
initmerge = pe.Node(fsl.Merge(dimension="t"),
                    name="initmerge")

# Mean the merged affine-transformed files
initmean = pe.Node(fsl.ImageMaths(op_string="-Tmean",
                                  suffix="_mean"),
                   name="initmean")

# Flip the mean image
initflip = pe.Node(fsl.SwapDimensions(new_dims=("-x","y","z")),
                   name="initflip")

# Add the mean images together
initadd = pe.Node(fsl.ImageMaths(op_string="-add"),
                  name="initadd")

# And then divide by 2 to preserve scale
inittemp = pe.Node(fsl.ImageMaths(op_string="-div 2"),
                  name="inittemp")

# Sink the initial template
initsink = pe.Node(nio.DataSink(base_directory=analysis_dir,
                                substitutions=[("mprage_gm_flirt_merged_mean_maths_maths",
                                                "initial_gm_template")]),
                   name="initsink")

# Define and connect the initial template flow
initialtemp = pe.Workflow(name="initialtemp", base_dir=working_dir)
flutil.archive_crashdumps(initialtemp)
initialtemp.connect([
    (initsource,    initmerge,   [("gray_image", "in_files")]),
    (initmerge,     initmean,    [("merged_file", "in_file")]),   
    (initmean,      initflip,    [("out_file", "in_file")]),
    (initmean,      initadd,     [("out_file", "in_file")]),
    (initflip,      initadd,     [("out_file", "in_file2")]),
    (initadd,       inittemp,    [("out_file", "in_file")]),
    (inittemp,      initsink,    [("out_file", "template.@initial_template")]),
    ])

# Set up a node to grab this initial template
initgrabber = pe.Node(nio.DataGrabber(outfields=["initial_template"],
                                      base_directory=analysis_dir,
                                      template="template/initial_gm_template.nii.gz"),
                      name="initgrabber")

# Now FLIRT each gray matter image to this initial template
flirttemp = pe.Node(fsl.FLIRT(searchr_x=[-180,180],
                                 searchr_y=[-180,180],
                                 searchr_z=[-180,180]),
                       name="flirttemp")

# Then improve the registration with FNIRT
fnirttemp = pe.Node(fsl.FNIRT(fieldcoeff_file=True,
                              config_file=pjoin(os.environ["FSLDIR"],"etc/flirtsch/GM_2_MNI152GM_2mm.cnf")),
                       name="fnirttemp")

# Apply the FNIRT warpfield with spline interpolation
warptemp = pe.Node(fsl.ApplyWarp(interp="spline"),
                   name="warptemp")

# Sink the nonlinear-transformed gray matter images
tempregsink = pe.Node(nio.DataSink(base_directory=analysis_dir,
                                   substitutions=[("_subject_id_", "")]),
                      name="tempregsink")

# Define and connect the main template registration flow
mainregflow = pe.Workflow(name="mainreg", base_dir=working_dir)
flutil.archive_crashdumps(mainregflow)
mainregflow.connect([
    (tempsubjects,  graysource,  [("subject_id", "subject_id")]),
    (initgrabber,   flirttemp,   [("initial_template", "reference")]),
    (graysource,    flirttemp,   [("gray_image", "in_file")]),
    (initgrabber,   fnirttemp,   [("initial_template", "ref_file")]),
    (flirttemp,     fnirttemp,   [("out_matrix_file", "affine_file")]),
    (graysource,    fnirttemp,   [("gray_image", "in_file")]),
    (initgrabber,   warptemp,    [("initial_template", "ref_file")]),
    (graysource,    warptemp,    [("gray_image", "in_file")]),
    (fnirttemp,     warptemp,    [("fieldcoeff_file", "field_file")]),
    (warptemp,      tempregsink, [("out_file", "@warped_gray")]),
    ])



# Grab the warped gray matter images
tempsource = pe.MapNode(nio.DataGrabber(infields=["subject_id"],
                                        outfields=["gray_image"],
                                        base_directory=analysis_dir,
                                        sort_filelist=True,
                                        template="%s/mprage_gm_warp.nii.gz"),
                        iterfield=["subject_id"],
                        name="tempsource")
tempsource.inputs.template_args=dict(gray_image=[["subject_id"]])
tempsource.inputs.subject_id = args.subjects


# Merge the nonlinear-registered images
tempmerge = pe.Node(fsl.Merge(dimension="t"),
                    name="tempmerge")

# Mean the merged nonlinear-transformed files
tempmean = pe.Node(fsl.ImageMaths(op_string="-Tmean",
                                  suffix="_mean"),
                   name="tempmean")

# Flip the mean image
tempflip = pe.Node(fsl.SwapDimensions(new_dims=("-x","y","z")),
                   name="tempflip")

# Add the mean images together
tempadd = pe.Node(fsl.ImageMaths(op_string="-add"),
                  name="tempadd")

# And then divide by 2 to preserve scale, producing the final template
maketemp = pe.Node(fsl.ImageMaths(op_string="-div 2"),
                  name="maketemp")

# Sink the template
tempsink = pe.Node(nio.DataSink(base_directory=analysis_dir,
                                substitutions=[("mprage_gm_warp_merged_mean_maths_maths",
                                                "final_gm_template")]),
                   name="tempsink")

# Define and connect the template workflow
maintemp = pe.Workflow(name="template", base_dir=working_dir)
flutil.archive_crashdumps(maintemp)
maintemp.connect([
    (tempsource,   tempmerge, [("gray_image", "in_files")]),
    (tempmerge,    tempmean,  [("merged_file", "in_file")]),
    (tempmean,     tempflip,  [("out_file", "in_file")]),
    (tempmean,     tempadd,   [("out_file", "in_file")]),
    (tempflip,     tempadd,   [("out_file", "in_file2")]),
    (tempadd,      maketemp,  [("out_file", "in_file")]),
    (maketemp,     tempsink,  [("out_file", "template.@template")]),
    ])
                                

# Normalization Workflow
# ----------------------

# Set up a node to grab the final template
tempgrabber = pe.Node(nio.DataGrabber(outfields=["final_template"],
                                      base_directory=analysis_dir,
                                      template="template/final_gm_template.nii.gz"),
                      name="initgrabber")

# Now FLIRT each gray matter image to this initial template
normflirt = pe.Node(fsl.FLIRT(searchr_x=[-180,180],
                                 searchr_y=[-180,180],
                                 searchr_z=[-180,180]),
                       name="normflirt")

# Then improve the registration with FNIRT
normfnirt = pe.Node(fsl.FNIRT(fieldcoeff_file=True,
                              jacobian_file=True,
                              config_file=pjoin(os.environ["FSLDIR"],"etc/flirtsch/GM_2_MNI152GM_2mm.cnf")),
                       name="normfnirt")

# Apply the FNIRT warpfield with spline interpolation
normwarp = pe.Node(fsl.ApplyWarp(interp="spline"),
                   name="normwarp")

# Modulate the warped file with its jacobian determinant
jacobmod = pe.Node(fsl.ImageMaths(op_string="-mul",
                                  suffix="_mod",
                                  out_data_type="float"),
                   name="jacobmod")

# Smooth the modulated images
kernel=3
smooth = pe.Node(fsl.ImageMaths(op_string="-s %.2f"%kernel,
                                suffix="_sm%d"%kernel),
                 name="smooth")

# Sink some images from this workflow
normsink = pe.Node(nio.DataSink(base_directory=analysis_dir,
                                substitutions=[("_subject_id_",""),
                                               ("warp", "norm")]),
                   name="normsink")

normalise = pe.Workflow(name="normalise", base_dir=working_dir)
flutil.archive_crashdumps(normalise)
normalise.connect([
    (subjsource,   graysource,  [("subject_id", "subject_id")]),
    (graysource,   normflirt,   [("gray_image", "in_file")]),
    (tempgrabber,  normflirt,   [("final_template", "reference")]),
    (normflirt,    normfnirt,   [("out_matrix_file", "affine_file")]),
    (graysource,   normfnirt,   [("gray_image", "in_file")]),
    (tempgrabber,  normfnirt,   [("final_template", "ref_file")]),
    (normfnirt,    normwarp,    [("fieldcoeff_file", "field_file")]),
    (graysource,   normwarp,    [("gray_image", "in_file")]),
    (tempgrabber,  normwarp,    [("final_template", "ref_file")]),
    (normwarp,     jacobmod,    [("out_file", "in_file")]),
    (normfnirt,    jacobmod,    [("jacobian_file", "in_file2")]),
    (jacobmod,     smooth,      [("out_file", "in_file")]),
    (normwarp,     normsink,    [("out_file", "@normalised_image")]),
    (jacobmod,     normsink,    [("out_file", "@modulated_image")]),
    (smooth,       normsink,    [("out_file", "@smoothed_image")]),
    ])

def workflow_runner(flow, stem):
    if any([a.startswith(stem) for a in args.workflows]) or args.workflows==["all"]:
        flow.run(inseries=args.inseries)

if __name__ == "__main__" and args.workflows is not None:
    # Run some things
    workflow_runner(seg, "seg")
    workflow_runner(initregflow, "temp")
    workflow_runner(initialtemp, "temp")
    workflow_runner(mainregflow, "temp")
    workflow_runner(maintemp, "temp")
    workflow_runner(normalise, "norm")
