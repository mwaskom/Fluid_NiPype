"""
Preprocessing for resting-state anaylsis
"""
import nibabel as nib
import nipype.interfaces.fsl as fsl       
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util   
import nipype.pipeline.engine as pe         
import nipype.algorithms.rapidart as ra      


preproc = pe.Workflow(name="preproc")

# Define the inputs for the preprocessing workflow
inputnode = pe.Node(util.IdentityInterface(fields=["timeseries",
                                                   "anatomical",
                                                   "warpfield",
                                                   "smooth_fwhm"]),
                    name="inputspec")

# Remove the first two frames to account for the lack of dummy acquisitions
trimmer = pe.MapNode(fsl.ExtractROI(t_min=2),
                     iterfield=["in_file"],
                     name="trimmer")

# Convert functional images to float representation
img2float = pe.MapNode(fsl.ImageMaths(out_data_type="float",
                                      op_string = "",
                                      suffix="_flt"),
                       iterfield=["in_file"],
                       name="img2float")

# Get the middle volume of each run for motion correction 
extractref = pe.MapNode(fsl.ExtractROI(t_size=1),
                         iterfield=["in_file"],
                         name = "extractref")

# Slice the example func for reporting
exampleslice = pe.MapNode(fsl.Slicer(image_width = 991,
                                     label_slices = False),
                          iterfield=["in_file"],
                          name="exampleslice")
exampleslice.inputs.sample_axial=2

# Motion correct to middle volume of each run
realign =  pe.MapNode(fsl.MCFLIRT(stages=4,
                                  interpolation="sinc",
                                  save_mats = True,
                                  save_plots = True,
                                  save_rms = True),
                      name="realign",
                      iterfield = ["in_file", "ref_file"])

# Plot the rotations, translations, and displacement parameters from MCFLIRT
plotrot = pe.MapNode(fsl.PlotMotionParams(in_source="fsl",
                                          plot_type="rotations"),
                     name="plotrotation", 
                     iterfield=["in_file"])

plottrans = pe.MapNode(fsl.PlotMotionParams(in_source="fsl",
                                            plot_type="translations"),
                       name="plottranslation", 
                       iterfield=["in_file"])

plotdisp = pe.MapNode(fsl.PlotMotionParams(in_source="fsl",
                                           plot_type="displacement"),
                      name="plotdisplacement",
                      iterfield=["in_file"])

# Perform slice-timing correction
slicetime = pe.MapNode(fsl.SliceTimer(interleaved=True,
                                      time_repetition=6),
                       iterfield=["in_file"],
                       name="slicetime")

# Get a mean image of the realigned timeseries
meanfunc1 = pe.MapNode(fsl.ImageMaths(op_string = "-Tmean",
                                      suffix="_mean"),
                       iterfield=["in_file"],
                       name="meanfunc1")

# Skullstrip the mean functional image
stripmean = pe.MapNode(fsl.BET(mask = True,
                               no_output=True,
                               frac = 0.3),
                       iterfield = ["in_file"],
                       name = "stripmean")

# Use the mask from skullstripping to strip each timeseries
maskfunc1 = pe.MapNode(fsl.ImageMaths(suffix="_bet",
                                      op_string="-mas"),
                       iterfield=["in_file", "in_file2"],
                       name = "maskfunc1")

# Define function for the first set of connections
def gettrimmedlength(func):
    """Return the desired length after removing two frames."""
    funcfile = func
    if isinstance(func, list):
        funcfile = func[0]
    _,_,_,timepoints = nib.load(funcfile).get_shape()
    return timepoints-2

def getmiddlevolume(func):
    """Return the middle volume index."""
    funcfile = func
    if isinstance(func, list):
        funcfile = func[0]
    _,_,_,timepoints = nib.load(funcfile).get_shape()
    return (timepoints/2)-1

# Connect the nodes for the first stage of preprocessing
preproc.connect([
    (inputnode,  trimmer,       [("timeseries", "in_file"),
                                 (("timeseries", gettrimmedlength), "t_size")]),
    (trimmer,    img2float,     [("roi_file", "in_file")]),
    (img2float,  extractref,    [("out_file", "in_file"), 
                                (("out_file", getmiddlevolume), "t_min")]),
    (extractref, exampleslice,  [("roi_file", "in_file")]),
    (img2float,  realign,       [("out_file", "in_file")]),
    (realign,  slicetime,       [("out_file", "in_file")]),
    (extractref, realign,       [("roi_file", "ref_file")]),
    (realign,    plotrot,       [("par_file", "in_file")]),
    (realign,    plottrans,     [("par_file", "in_file")]),
    (realign,    plotdisp,      [("rms_files","in_file")]),
    (slicetime,  meanfunc1,     [("slice_time_corrected_file", "in_file")]),
    (meanfunc1,  stripmean,     [("out_file", "in_file")]),
    (slicetime,  maskfunc1,     [("slice_time_corrected_file", "in_file")]),
    (stripmean,  maskfunc1,     [("mask_file", "in_file2")]),
    ])

# Determine the 2nd and 98th percentile intensities of each run
getthresh = pe.MapNode(fsl.ImageStats(op_string="-p 2 -p 98"),
                       iterfield = ["in_file"],
                       name="getthreshold")

# Threshold the first fun of the functional data at 10% of the 98th percentile
threshold = pe.MapNode(fsl.ImageMaths(out_data_type="char",
                                                suffix="_thresh"),
                       iterfield = ["in_file"],
                       name="threshold")


# Determine the median value of the functional runs using the mask
medianval = pe.MapNode(fsl.ImageStats(op_string="-k %s -p 50"),
                       iterfield = ["in_file", "mask_file"],
                       name="medianval")

# Dilate the mask
dilatemask = pe.MapNode(fsl.ImageMaths(suffix="_dil",
                                                 op_string="-dilF"),
                        iterfield=["in_file"],
                        name="dilatemask")

# Mask the runs again with this new mask
maskfunc2 = pe.MapNode(fsl.ImageMaths(suffix="_mask",
                                                op_string="-mas"),
                      iterfield=["in_file", "in_file2"],
                      name="maskfunc2")

# Get a new mean image from each functional run
meanfunc2 = pe.MapNode(fsl.ImageMaths(op_string="-Tmean",
                                                suffix="_mean"),
                       iterfield=["in_file"],
                       name="meanfunc2")

# Slice the mean func for reporting
meanslice = pe.MapNode(fsl.Slicer(image_width = 991,
                                  label_slices = False),
                       iterfield=["in_file"],
                       name="meanslice")
meanslice.inputs.sample_axial = 2


# Define a function for the second set of connections
def getthreshop(thresh):
    """Return an fslmaths op string to get10% of the intensity"""
    return "-thr %.10f -Tmin -bin"%(0.1*thresh[0][1])

# Connect the pipeline up through the new masked runs
preproc.connect([
    (maskfunc1,  getthresh,  [("out_file", "in_file")]),
    (getthresh,  threshold,  [(("out_stat", getthreshop), "op_string")]),
    (maskfunc1,  threshold,  [("out_file", "in_file")]),
    (slicetime,  medianval,  [("slice_time_corrected_file", "in_file")]),
    (threshold,  medianval,  [("out_file", "mask_file")]),
    (threshold,  dilatemask, [("out_file", "in_file")]),
    (slicetime,  maskfunc2,  [("slice_time_corrected_file", "in_file")]),
    (dilatemask, maskfunc2,  [("out_file", "in_file2")]),
    (maskfunc2,  meanfunc2,  [("out_file", "in_file")]),
    (meanfunc2,  meanslice,  [("out_file", "in_file")]),
    ])

# Use RapidART to detect motion/intensity outliers
art = pe.MapNode(ra.ArtifactDetect(use_differences = [True, False],
                                             use_norm = True,
                                             zintensity_threshold = 3,
                                             norm_threshold = 1,
                                             parameter_source = "FSL",
                                             mask_type = "file"),
                 iterfield=["realignment_parameters","realigned_files","mask_file"],
                 name="art")

# Plot a timecourse of the global mean intensity
plotmean = pe.MapNode(fsl.PlotTimeSeries(title="Global Mean Intensity"),
                      iterfield=["in_file"],
                      name="plotmean")

# Make connections to ART
preproc.connect([
    (realign,    art, [("par_file", "realignment_parameters")]),
    (maskfunc2,  art, [("out_file", "realigned_files")]),
    (dilatemask, art, [("out_file", "mask_file")]),
    (art,        plotmean, [("intensity_files", "in_file")]),   
    ])

# Merge the median values with the mean functional images into a coupled list
mergenode = pe.Node(util.Merge(2, axis="hstack"),
                    name="merge")

# Smooth in the volume with SUSAN
smooth = pe.MapNode(fsl.SUSAN(),
                    iterfield=["in_file", "brightness_threshold", "usans"],
                    name="smooth")

# Mask the smoothed data with the dilated mask
masksmoothfunc = pe.MapNode(fsl.ImageMaths(suffix="_mask",
                                                     op_string="-mas"),
                            iterfield=["in_file","in_file2"],
                            name="masksmoothfunc")
# SUSAN functions
def getbtthresh(medianvals):
    """Get the brightness threshold for SUSAN."""
    return [0.75*val for val in medianvals]

def getusans(inlist):
    """Return the usans at the right threshold."""
    return [[tuple([val[0],0.75*val[1]])] for val in inlist]

# Make the smoothing connections
preproc.connect([
    (meanfunc2,  mergenode,      [("out_file", "in1")]),
    (medianval,  mergenode,      [("out_stat", "in2")]),
    (maskfunc2,  smooth,         [("out_file", "in_file")]),
    (inputnode,  smooth,         [("smooth_fwhm", "fwhm")]),
    (medianval,  smooth,         [(("out_stat", getbtthresh), "brightness_threshold")]),
    (mergenode,  smooth,         [(("out", getusans), "usans")]),
    (smooth,     masksmoothfunc, [("smoothed_file", "in_file")]),
    (dilatemask, masksmoothfunc, [("out_file", "in_file2")]),
    ])

# Unzip the anatomical so it plays nicely with Conn
unzip_anat = pe.Node(fs.MRIConvert(out_type="nii"),
                     name="unzipanatomical")

# Make unzipping connections
preproc.connect(inputnode, "anatomical", unzip_anat, "in_file")

# Define the outputs of the preprocessing that will be used by the model
outputnode = pe.Node(util.IdentityInterface(fields=["smoothed_timeseries",
                                                    "unsmoothed_timeseries",
                                                    "example_func",
                                                    "mean_func",
                                                    "T1_warped",
                                                    "functional_mask",
                                                    "realignment_parameters",
                                                    "outlier_files"]),
                     name="outputspec")

# Define the outputs that will go into a report
reportnode = pe.Node(util.IdentityInterface(fields=["example_func",
                                                    "mean_func",
                                                    "rotation_plot",
                                                    "outlier_volumes",
                                                    "intensity_plot",
                                                    "translation_plot",
                                                    "displacement_plot"]),
                     name="report")

# Make connections to the output nodes
preproc.connect([
    (maskfunc2,      outputnode, [("out_file", "unsmoothed_timeseries")]),
    (masksmoothfunc, outputnode, [("out_file", "smoothed_timeseries")]),
    (unzip_anat,     outputnode, [("out_file", "T1_warped")]),
    (extractref,     outputnode, [("roi_file", "example_func")]),
    (meanfunc2,      outputnode, [("out_file", "mean_func")]),
    (dilatemask,     outputnode, [("out_file", "functional_mask")]),
    (realign,        outputnode, [("par_file", "realignment_parameters")]),
    (art,            outputnode, [("outlier_files", "outlier_files")]),
    (art,            reportnode, [("outlier_files", "outlier_volumes")]),
    (exampleslice,   reportnode, [("out_file", "example_func")]),
    (plotmean,       reportnode, [("out_file", "intensity_plot")]),    
    (meanslice,      reportnode, [("out_file", "mean_func")]),
    (plotrot,        reportnode, [("out_file", "rotation_plot")]),
    (plottrans,      reportnode, [("out_file", "translation_plot")]),
    (plotdisp,       reportnode, [("out_file", "displacement_plot")]),
    ])
