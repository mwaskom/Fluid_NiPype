#! /usr/bin/env python
"""
    Use FLIRT and FNIRT to register images from recon-all to the FSL target.
"""

import os
import sys
import nipype.interfaces.fsl as fsl
import nipype.interfaces.freesurfer as fs
from nipype.utils.filemanip import FileNotFoundError

if len(sys.argv) < 2:
    sys.exit("USAGE: fluid_register.py [-f] SUBJECT")

subject = sys.argv[-1]
force = False
if "-f" in sys.argv:
    force = True

# Get target images
target_brain = fsl.Info.standard_image("avg152T1_brain.nii.gz")
target_head =  fsl.Info.standard_image("avg152T1.nii.gz")
fnirt_cfg = "/usr/share/fsl/4.1/etc/flirtsch/T1_2_MNI152_2mm.cnf"

# Get top-level dir
subjdir = os.path.join("/mindhive/gablab/fluid/Data", subject)

# Get recon images
brainmask_mgz = os.path.join(subjdir, "mri/brainmask.mgz")
t1_mgz = os.path.join(subjdir, "mri/T1.mgz")

# Define the images to write
brainmask = "brainmask.nii.gz"
t1 = "T1.nii.gz"

brainmask_flirted = "brainmask_flirted.nii.gz"
t1_fnirted = "T1_fnirted.nii.gz"

flirtmat = "affine.mat"
fnirtfield = "warpfield.nii.gz"

qcpng = "final_reg.png"

# Set up the output dir
regdir = os.path.join(subjdir, "registration")
try:
    os.mkdir(regdir)
except:
    pass

origdir = os.getcwd()
os.chdir(regdir)

# Set up logging
logfile = open("registration.log","w")

def log(interface, result):
    logfile.write("\n".join([interface.cmdline,
                             result.runtime.stdout,
                             result.runtime.stderr]))

try:
    # Brainmask to nifti
    cvt = fs.MRIConvert(in_file=brainmask_mgz, out_file=brainmask)
    if force or not os.path.exists(brainmask):
        res = cvt.run()
        log(cvt, res)

    # T1 to nifti
    cvt = fs.MRIConvert(in_file=t1_mgz, out_file=t1)
    if force or not os.path.exists(t1):
        res = cvt.run()
        log(cvt, res)

    # FLIRT brainmask to mni152_brain
    flirt = fsl.FLIRT(in_file         = brainmask, 
                      reference       = target_brain,
                      out_file        = brainmask_flirted,
                      out_matrix_file = flirtmat,
                      searchr_x       = [-180, 180],
                      searchr_y       = [-180, 180],
                      searchr_z       = [-180, 180])
    if force or (not os.path.exists(brainmask_flirted) or not os.path.exists(flirtmat)):
        res = flirt.run()
        log(flirt, res)

    # FNIRT T1 to mni152
    fnirt = fsl.FNIRT(in_file         = t1,
                      ref_file        = target_head,
                      config_file     = fnirt_cfg,
                      affine_file     = flirtmat,
                      warped_file     = t1_fnirted,
                      fieldcoeff_file = fnirtfield)
    if force or (not os.path.exists(t1_fnirted) or not os.path.exists(fnirtfield)):
        try:
            res = fnirt.run()
        # FNIRT is currently being weird about writing
        # files with names other than the defaults that
        # get generated, so we'll deal with that here
        except FileNotFoundError:
            if not os.path.exists(fnirtfield):
                raise FileNotFoundError("%s was not written by FNIRT")
        log(fnirt, res)

    # Move the FNIRT log contents into the main log
    fnirtlog = "T1_to_avg152T1.log"
    if os.path.exists(fnirtlog):
        logfile.write(open(fnirtlog).read())
        os.remove(fnirtlog)

    # Slice output for qc
    slicer = fsl.Slicer(in_file      = t1_fnirted,
                        image_edges  = target_brain,
                        image_width  = 750,
                        sample_axial = 3,
                        out_file     = qcpng)
    if force or not os.path.exists(qcpng):
        res = slicer.run()
        log(slicer, res)

except:
    # Deal with exceptions
    logfile.close()
    os.chdir(origdir)
    raise
