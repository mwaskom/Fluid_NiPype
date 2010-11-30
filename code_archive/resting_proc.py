"""
    Preprocessing module for Fluid Intelligence resting state
"""


import os                               
from glob import glob
import nipype.interfaces.io as nio       
import nipype.interfaces.fsl as fsl       
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util   
import nipype.pipeline.engine as pe         
import nipype.algorithms.rapidart as ra      
import nibabel as nib

preproc = pe.Workflow(name="preproc",
    base_dir="/mindhive/gablab/fluid/Analysis/NiPype/workingdir/resting")

infosource = pe.Node(interface=util.IdentityInterface(fields=["subject_id"]),
                     name="infosource")

data_dir = "/mindhive/gablab/fluid/Data"
subject_list = [p.split("/")[-3] for p in glob(os.path.join(data_dir,"gf??/bold/Resting.nii.gz"))]
subject_list = [s for s in subject_list if 
                    os.path.exists(os.path.join(data_dir,s,"registration","warpfield.nii.gz"))]
print "Subjects:", " ".join(subject_list)

infosource.iterables = ("subject_id", subject_list)

datasource = pe.Node(interface=nio.DataGrabber(infields=["subject_id"],
                                               outfields=["resting_timeseries","warpfield"]),
                     name="datasource")

datasource.inputs.base_directory = data_dir
datasource.inputs.template = "%s/%s/%s.nii.gz"
datasource.inputs.template_args = dict(resting_timeseries=[["subject_id","bold","Resting"]],
                                       warpfield=[["subject_id","registration","warpfield"]])

# Define the inputs for the preprocessing workflow
inputnode = pe.Node(interface=util.IdentityInterface(fields=["resting_timeseries",
                                                             "warpfield",
                                                             "subject_id"]),
                    name="inputspec")

preproc.connect([
    (infosource, datasource, [("subject_id", "subject_id")]),
    (infosource, inputnode,  [("subject_id", "subject_id")]),
    (datasource, inputnode,  [("resting_timeseries", "resting_timeseries")]),
    (datasource, inputnode,  [("warpfield", "warpfield")]),
    ])

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

slicetime = pe.MapNode(interface=fsl.SliceTimer(interleaved=True,
                                                time_repetition=6),
                       iterfield=["in_file"],
                       name="slicetime")

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
    (inputnode,  img2float,   [("resting_timeseries", "in_file")]),
    (img2float,  extractref,  [("out_file", "in_file"), 
                              (("out_file", getmiddlevolume), "t_min")]),
    (img2float,  slicetime,   [("out_file", "in_file")]),
    (slicetime,  realign,     [("slice_time_corrected_file", "in_file")]),
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
# Connect the final steps of the volume preprocessing
preproc.connect([
    (volsmooth,    volmaskfunc,  [("smoothed_file", "in_file")]),
    (dilatemask,   volmaskfunc,  [("out_file", "in_file2")]),
    ])

# Register each mean functional to Freesurfer anat space
regfunc = pe.MapNode(fs.BBRegister(contrast_type="t2",
                                   init="fsl",
                                   out_fsl_file=True),
                     iterfield=["source_file"],
                     name="regfunc")

# Apply the nonlinear warp to the timeseries
warpfunc = pe.MapNode(fsl.ApplyWarp(ref_file=fsl.Info.standard_image("avg152T1.nii.gz")),
                      iterfield=["in_file","premat"],
                      name="warpfunc")

# Mean the warped timeseries for reporting
meanwarp = pe.MapNode(fsl.ImageMaths(op_string="-Tmean",
                                     suffix="_mean"),
                      iterfield=["in_file"],
                      name="meanwarp")

# Slice the mean warped image 
slicewarp = pe.MapNode(fsl.Slicer(middle_slices=True,
                                  label_slices=False,
                                  image_edges=fsl.Info.standard_image("avg152T1_brain.nii.gz")),
                       iterfield=["in_file"],
                       name="slicewarp")

# Unzip the final normalized timeseries for SPM
unzip = pe.MapNode(fs.MRIConvert(out_type="nii"),
                   iterfield=["in_file"],
                   name="unzip")

# Connect the registration pipeline
preproc.connect([
    (inputnode,   regfunc,   [("subject_id", "subject_id")]),
    (meanfunc2,   regfunc,   [("out_file", "source_file")]),
    (volmaskfunc, warpfunc,  [("out_file", "in_file")]),
    (inputnode,   warpfunc,  [("warpfield", "field_file")]),
    (regfunc,     warpfunc,  [("out_fsl_file", "premat")]),
    (warpfunc,    meanwarp,  [("out_file", "in_file")]),
    (warpfunc,    unzip,     [("out_file", "in_file")]),
    (meanwarp,    slicewarp, [("out_file", "in_file")]),
    ]) 

sinkrest = pe.Node(interface=nio.DataSink(), name="sinkrest")
sinkrest.inputs.base_directory="/mindhive/gablab/fluid/Analysis/NiPype/resting"
sinkrest.inputs.parameterization=False

preproc.connect([
    (infosource, sinkrest, [("subject_id", "container"),
                            (("subject_id", lambda x: "_subject_id_%s"%x), "strip_dir")]),
    (unzip,      sinkrest, [("out_file", "@resting_ts")]),
    (slicewarp,  sinkrest, [("out_file", "@warppng")]),
    ])
