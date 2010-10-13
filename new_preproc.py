"""
    NEW - Normalization preproc
    Preprocessing module for Fluid Intelligence fMRI paradigms.
"""


import os                               
import sys
import nipype.interfaces.io as nio       
import nipype.interfaces.fsl as fsl       
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.spm as spm
import nipype.interfaces.utility as util   
import nipype.pipeline.engine as pe         
import nipype.algorithms.rapidart as ra      
import nipype.algorithms.modelgen as model
import nibabel as nib

preproc = pe.Workflow(name="preproc")

# Define the inputs for the preprocessing workflow
inputnode = pe.Node(interface=util.IdentityInterface(fields=["func",
                                                             "subject_id"]),
                    name="inputspec")

# Convert functional images to float representation
img2float = pe.MapNode(interface=fsl.ImageMaths(out_data_type="float",
                                             op_string = "",
                                             suffix="_dtype"),
                       iterfield=["in_file"],
                       name="img2float")

# Get the middle volume of each run for motion correction 
extractref = pe.MapNode(interface=fsl.ExtractROI(t_size=1),
                         iterfield=["in_file"],
                         name = "extractref")

# Motion correct to middle volume of each run
realign =  pe.MapNode(interface=fsl.MCFLIRT(save_mats = True,
                                            save_plots = True,
                                            save_rms = True),
                      name="realign",
                      iterfield = ["in_file", "ref_file"])

# Plot the rotations and translation parameters from MCFLIRT
plotmotion = pe.MapNode(interface=fsl.PlotMotionParams(in_source="fsl"),
                        name="plotmotion", 
                        iterfield=["in_file"],
                        iterables = ("plot_type", ["rotations","translations"]))

# Plot the mean displacement parameters from MCFLIRT
plotdisp = pe.MapNode(interface=fsl.PlotMotionParams(in_source="fsl",
                                                     plot_type="displacement"),
                      name="plotdisplacement",
                      iterfield=["in_file"])

# Get a mean image of the realigned timeseries
meanfunc1 = pe.MapNode(interface=fsl.ImageMaths(op_string = "-Tmean",
                                               suffix="_mean"),
                       iterfield=["in_file"],
                       name="meanfunc1")

# Skullstrip the mean functional image
stripmean = pe.MapNode(interface=fsl.BET(mask = True,
                                         no_output=True,
                                         frac = 0.3),
                       iterfield = ["in_file"],
                       name = "stripmean")

# Use the mask from skullstripping to strip each timeseries
maskfunc1 = pe.MapNode(interface=fsl.ImageMaths(suffix="_bet",
                                               op_string="-mas"),
                       iterfield=["in_file", "in_file2"],
                       name = "maskfunc1")

# Define function for the first set of connections
def getmiddlevolume(func):
    """Return the middle volume index."""
    funcfile = func
    if isinstance(func, list):
        funcfile = func[0]
    _,_,_,timepoints = nib.load(funcfile).get_shape()
    return (timepoints/2)-1

# Connect the nodes for the first stage of preprocessing
preproc.connect([
    (inputnode,  img2float,   [("func", "in_file")]),
    (img2float,  extractref,  [("out_file", "in_file"), 
                              (("out_file", getmiddlevolume), "t_min")]),
    (img2float,  realign,     [("out_file", "in_file")]),
    (extractref, realign,     [("roi_file", "ref_file")]),
    (realign,    plotmotion,  [("par_file", "in_file")]),
    (realign,    plotdisp,    [("rms_files","in_file")]),
    (realign,    meanfunc1,   [("out_file", "in_file")]),
    (meanfunc1,  stripmean,   [("out_file", "in_file")]),
    (realign,    maskfunc1,   [("out_file", "in_file")]),
    (stripmean, maskfunc1,    [("mask_file", "in_file2")]),
    ])

# Determine the 2nd and 98th percentile intensities of each run
getthresh = pe.MapNode(interface=fsl.ImageStats(op_string="-p 2 -p 98"),
                       iterfield = ["in_file"],
                       name="getthreshold")

# Threshold the first fun of the functional data at 10% of the 98th percentile
threshold = pe.MapNode(interface=fsl.ImageMaths(out_data_type="char",
                                                suffix="_thresh"),
                       iterfield = ["in_file"],
                       name="threshold")


# Determine the median value of the functional runs using the mask
medianval = pe.MapNode(interface=fsl.ImageStats(op_string="-k %s -p 50"),
                       iterfield = ["in_file", "mask_file"],
                       name="medianval")

# Dilate the mask
dilatemask = pe.MapNode(interface=fsl.ImageMaths(suffix="_dil",
                                                 op_string="-dilF"),
                        iterfield=["in_file"],
                        name="dilatemask")

# Mask the runs again with this new mask
maskfunc2 = pe.MapNode(interface=fsl.ImageMaths(suffix="_mask",
                                                op_string="-mas"),
                      iterfield=["in_file", "in_file2"],
                      name="maskfunc2")

# Get a new mean image from each functional run
meanfunc2 = pe.MapNode(interface=fsl.ImageMaths(op_string="-Tmean",
                                                suffix="_mean"),
                       iterfield=["in_file"],
                       name="meanfunc2")

# Merge the median values with the mean functional images into a coupled list
mergenode = pe.Node(interface=util.Merge(2, axis="hstack"),
                    name="merge")

# Define a function for the second set of connections
def getthreshop(thresh):
    """Return an fslmaths op string to get10% of the intensity"""
    return "-thr %.10f -Tmin -bin"%(0.1*thresh[0][1])

# Connect the pipeline up through the new masked runs
preproc.connect([
    (maskfunc1,  getthresh,  [("out_file", "in_file")]),
    (getthresh,  threshold,  [(("out_stat", getthreshop), "op_string")]),
    (maskfunc1,  threshold,  [("out_file", "in_file")]),
    (realign,    medianval,  [("out_file", "in_file")]),
    (threshold,  medianval,  [("out_file", "mask_file")]),
    (threshold,  dilatemask, [("out_file", "in_file")]),
    (realign,    maskfunc2,  [("out_file", "in_file")]),
    (dilatemask, maskfunc2,  [("out_file", "in_file2")]),
    (maskfunc2,  meanfunc2,  [("out_file", "in_file")]),
    (meanfunc2,  mergenode,  [("out_file", "in1")]),
    (medianval,  mergenode,  [("out_stat", "in2")]),
    ])

# Use RapidART to detect motion/intensity outliers
art = pe.MapNode(interface=ra.ArtifactDetect(use_differences = [True, False],
                                             use_norm = True,
                                             zintensity_threshold = 3,
                                             norm_threshold = 1,
                                             parameter_source = "FSL",
                                             mask_type = "file"),
                 iterfield=["realignment_parameters","realigned_files","mask_file"],
                 name="art")

# Make connections to ART
preproc.connect([
    (realign,    art, [("par_file", "realignment_parameters")]),
    (maskfunc2,  art, [("out_file", "realigned_files")]),
    (dilatemask, art, [("out_file", "mask_file")]),
    ])

# Set volume smoothing kernel size
volsmoothval = pe.Node(util.IdentityInterface(fields=["fwhm"]),
                       iterables = ("fwhm", [5.]),
                       name="volsmoothval")


# Smooth in the volume with SUSAN
volsmooth = pe.MapNode(interface=fsl.SUSAN(),
                       iterfield=["in_file", "brightness_threshold", "usans"],
                       name="volsmooth")

# SUSAN functions
def getbtthresh(medianvals):
    """Get the brightness threshold for SUSAN."""
    return [0.75*val for val in medianvals]

def getusans(inlist):
    """Return the usans at the right threshold."""
    return [[tuple([val[0],0.75*val[1]])] for val in inlist]

# Make connections to SUSAN
preproc.connect([
    (maskfunc2,    volsmooth, [("out_file", "in_file")]),
    (volsmoothval, volsmooth, [("fwhm", "fwhm")]),
    (medianval,    volsmooth, [(("out_stat", getbtthresh), "brightness_threshold")]),
    (mergenode,    volsmooth, [(("out", getusans), "usans")]),
    ])

# Mask the smoothed data with the dilated mask
volmaskfunc = pe.MapNode(interface=fsl.ImageMaths(suffix="_mask",
                                                  op_string="-mas"),
                         iterfield=["in_file","in_file2"],
                         name="volmaskfunc")

# Scale the median value of the run to 10000
volmeanscale = pe.MapNode(interface=fsl.ImageMaths(suffix="_gms"),
                          iterfield=["in_file","op_string"],
                          name="volmeanscale")

# High-pass filter the timeseries
volhighpass = pe.MapNode(interface=fsl.ImageMaths(suffix="_vol_hpf",
                                                  op_string="-bptf 64 -1"),
                         iterfield=["in_file"],
                         name="volhighpass")

# Define a function for intensity normalization
def getmeanscaleop(medianvals):
    """Get an fslmaths op string for intensity normalization."""
    return ["-mul %.10f"%(10000./val) for val in medianvals]

# Connect the final steps of the volume preprocessing
preproc.connect([
    (volsmooth,    volmaskfunc,  [("smoothed_file", "in_file")]),
    (dilatemask,   volmaskfunc,  [("out_file", "in_file2")]),
    (volmaskfunc,  volmeanscale, [("out_file", "in_file")]),
    (medianval,    volmeanscale, [(("out_stat", getmeanscaleop), "op_string")]),
    (volmeanscale, volhighpass,  [("out_file", "in_file")]),
    ])

# Register each mean functional to Freesurfer anat space
regfunc = pe.MapNode(fs.BBRegister(contrast_type="t2",
                                   init="fsl",
                                   out_fsl_file=True),
                     iterfield=["source_file"],
                     name="regfunc")

# XXX Insert slicer report node here

# Get the T1 and brainmask from recon-all
fssource = pe.Node(nio.FreeSurferSource(subjects_dir="/mindhive/gablab/fluid/Data"),
                   name="fssource")

# Convert the brainmask to Nifti so FLIRT can read it
niftimask = pe.Node(fs.MRIConvert(out_type="niigz"),
                    name="niftimask")

# Get the standard-space FLIRT target
targetbrain = fsl.Info.standard_image("avg152T1_brain.nii.gz")

# Register the brainmask to the target using a 12-dof affine transformation
regstruct = pe.Node(fsl.FLIRT(reference=targetbrain,
                              searchr_x=[-180,180],
                              searchr_y=[-180,180],
                              searchr_z=[-180,180]),
                    name="regstruct")

# XXX Insert slicer report node here
# Convert the T1 to nifti so FNIRT can read it
niftit1 = pe.Node(fs.MRIConvert(out_type="niigz"),
                  name="niftit1")

# Path to the standard FNIRT config file
fnirtcfg = "/usr/share/fsl/4.1/etc/flirtsch/T1_2_MNI152_2mm.cnf"

# Get the standard space fnirt target
targethead = fsl.Info.standard_image("avg152T1.nii.gz")

# Use FNIRT to get a nonlinear transformation to the target
fnirt = pe.Node(fsl.FNIRT(config_file=fnirtcfg,
                          ref_file=targethead),
                          field_file=True,
                name="fnirt")

# Apply the nonlinear warp to the timeseries
warpfunc = pe.MapNode(fsl.ApplyWarp(ref_file=targethead),
                      iterfield=["in_file","premat"],
                      name="warpfunc")

# Concatenate the func to anat and anat to standard transform matrices
matconcat = pe.MapNode(fsl.ConvertXFM(concat_xfm=True),
                       iterfield=["in_file"],
                       name="matconcat")

# Apply the concatenated transformation to each timeseries
funcxfm = pe.MapNode(fsl.FLIRT(apply_xfm=True,
                               reference=targetbrain),
                     iterfield=["in_file", "in_matrix_file"],
                     name="funcxfm")

# Connect the registration pipeline
preproc.connect([
    (inputnode,   fssource,  [("subject_id", "subject_id")]),
    (inputnode,   regfunc,   [("subject_id", "subject_id")]),
    (meanfunc2,   regfunc,   [("out_file", "source_file")]),
    (fssource,    niftimask, [("brainmask", "in_file")]),
    (niftimask,   regstruct, [("out_file", "in_file")]),
    (fssource,    niftit1,   [("T1", "in_file")]),
    (niftit1,     fnirt,     [("out_file", "in_file")]),
    (regstruct,   fnirt,     [("out_matrix_file", "affine_file")]),
    (regfunc,     matconcat, [("out_fsl_file", "in_file")]),
    (regstruct,   matconcat, [("out_matrix_file", "in_file2")]),
    (fnirt,       warpfunc,  [("fieldcoeff_file", "field_file")]),
    (regfunc,     warpfunc,  [("out_fsl_file", "premat")]),
    (volhighpass, warpfunc,  [("out_file", "in_file")]),
    (matconcat,   funcxfm,   [("out_file", "in_matrix_file")]),
    (volhighpass, funcxfm,   [("out_file", "in_file")]),
    ])
