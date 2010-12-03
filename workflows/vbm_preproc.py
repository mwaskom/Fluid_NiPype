"""
A workflow for generic VBM-style preprocessing of statistical images.
This was written to process FA maps generated from DTI scans, but it
should work with arbitrary statistical images.

The steps are:
- warping to MNI space
- masking with an MNI template brain mask
- smoothing with an isotropic kernel

Input node takes the following inputs:
- A statistical image (the image to process)
- FSL style affine matrix to coregister the statistical image
  to an anatomical scan.
- A nonlinear warpfield that gives the transformation between 
  the native anatomical and the MNI brain.
- Smoothing kernel size in mm FWHM

Output node returns a processed image
Report node returns pngs for quality control on the registration

"""
from __future__ import division 

import numpy as np
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe

preproc = pe.Workflow(name="preproc")

# Define the inputs to the workflow
inputnode = pe.Node(util.IdentityInterface(fields=["stat_image",
                                                   "affine_mat",
                                                   "warpfield",
                                                   "smooth_fwhm"]),
                    name="inputspec")

# Warp the stat image to FSL's standard MNI template (with 2x2x2 voxels)
mni_brain = fsl.Info.standard_image("avg152T1_brain.nii.gz") 

applywarp = pe.Node(fsl.ApplyWarp(ref_file=mni_brain),
                    name="applywarp")

# Slice the warp registration for quality control
slicewarp = pe.Node(fsl.Slicer(label_slices=False,
                               show_orientation=False,
                               image_edges=mni_brain,
                               middle_slices=True),
                    name="slicewarp")
# Mask the FA image with the template mask
mni_mask = fsl.Info.standard_image('MNI152_T1_2mm_brain_mask_dil.nii.gz')

applymask = pe.Node(fsl.ImageMaths(suffix="_mask",
                                   op_string="-mas",
                                   in_file2=mni_mask),
                    name="applymask")


# Isotropically smooth the masked image
smooth = pe.Node(fsl.ImageMaths(suffix="_smooth"),
                 name="smooth")

# Slice the final image for quick quality control
sliceimg = pe.Node(fsl.Slicer(image_width=910,
                              label_slices=False),
                   name="sliceimg")
sliceimg.inputs.sample_axial=2

outputnode = pe.Node(util.IdentityInterface(fields=["stat_image"]),
                     name="outputspec")

report = pe.Node(util.IdentityInterface(fields=["reg_slices", "image_slices"]),
                 name="report")


# Define a function to convert a fwhm value into a fslmaths smoothing op_string
def getsmoothop(fwhm):
    """Convert a kernel fwhm to sigma and return an fslmaths op string"""
    sigma = fwhm/np.sqrt(8 * np.log(2)) 
    return "-s %.8f"%sigma

preproc.connect([
    (inputnode, applywarp,  [("stat_image", "in_file"),
                             ("warpfield", "field_file"),
                             ("affine_mat", "premat")]),
    (inputnode, smooth,     [(("smooth_fwhm", getsmoothop), "op_string")]),
    (applywarp, slicewarp,  [("out_file", "in_file")]),
    (applywarp, applymask,  [("out_file", "in_file")]),
    (applymask, smooth,     [("out_file", "in_file")]),
    (smooth,    sliceimg,   [("out_file", "in_file")]),
    (smooth,    outputnode, [("out_file", "stat_image")]),
    (slicewarp, report,     [("out_file", "reg_slices")]),
    (sliceimg,  report,     [("out_file", "image_slices")]),
    ])
