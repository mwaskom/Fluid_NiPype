"""
Preprocessing module for Fluid Intelligence fMRI paradigms.

Input spec node takes three inputs:
    - Timeseries (image files)
    - Voxel shift map
    - Fieldmap brain (forward-warped skull-stripped magnitude volume)
    - Highpass filter cutoff (in TRs)
    - FWHM of smoothing kernel for SUSAN (in mms)

Output node returns these files:
    - Smoothed timeseries (fully preprocessed and smoothed timeseries in native space)
    - Unsmoothed timeseries (identical steps except no smoothing in the volume)
    - Example func (target volume for MCFLIRT realignment)
    - Mean func (unsmoothed mean functional)
    - Funcational mask (binary dilated brainmask in functional space)
    - Realignment parameters (text files from MCFLIRT)
    - Outlier Files (outlier text files from ART)

Reporting node returns these files
    - Plotted estimated rotations from MCFLIRT
    - Plotted estimated translastion from MCFLIRT
    - Plotted estimated relative and absolute displacement from MCFLIRT
    - Plotted global mean intensity value
    - Sliced png of the example func (MCFLIRT target)
    - Sliced png of the unsmoothed mean functional volume

"""
import nibabel as nib
import nipype.interfaces.fsl as fsl       
import nipype.interfaces.utility as util   
import nipype.pipeline.engine as pe         
import nipype.algorithms.rapidart as ra      

preproc = pe.Workflow(name="preproc")

# Define the inputs for the preprocessing workflow
inputnode = pe.Node(interface=util.IdentityInterface(fields=["timeseries",
                                                             "voxel_shift_map",
                                                             "fieldmap_brain",
                                                             "hpf_cutoff",
                                                             "smooth_fwhm"]),
                    name="inputspec")

# Convert functional images to float representation
img2float = pe.MapNode(interface=fsl.ImageMaths(out_data_type="float",
                                             op_string = "",
                                             suffix="_flt"),
                       iterfield=["in_file"],
                       name="img2float")

# Get the middle volume of each run for motion correction 
extractref = pe.MapNode(interface=fsl.ExtractROI(t_size=1),
                         iterfield=["in_file"],
                         name = "extractref")

# Slice the example func for reporting
exampleslice = pe.MapNode(fsl.Slicer(image_width = 572,
                                     label_slices = False),
                          iterfield=["in_file"],
                          name="exampleslice")
exampleslice.inputs.sample_axial=2

# Motion correct to middle volume of each run
realign =  pe.MapNode(interface=fsl.MCFLIRT(stages=4,
					    interpolation="sinc",
					    save_mats = True,
                                            save_plots = True,
                                            save_rms = True),
                      name="realign",
                      iterfield = ["in_file", "ref_file"])

# Plot the rotations, translations, and displacement parameters from MCFLIRT
plotrot = pe.MapNode(interface=fsl.PlotMotionParams(in_source="fsl",
                                                    plot_type="rotations"),
                     name="plotrotation", 
                     iterfield=["in_file"])

plottrans = pe.MapNode(interface=fsl.PlotMotionParams(in_source="fsl",
                                                    plot_type="translations"),
                       name="plottranslation", 
                       iterfield=["in_file"])

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

# Use the mask to strip the example func
maskexample = pe.MapNode(fsl.ImageMaths(suffix="_bet",
                                        op_string="-mas"),
                         iterfield=["in_file", "in_file2"],
                         name = "maskexample")

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
    (inputnode,  img2float,     [("timeseries", "in_file")]),
    (img2float,  extractref,    [("out_file", "in_file"), 
                                (("out_file", getmiddlevolume), "t_min")]),
    (extractref, exampleslice,  [("roi_file", "in_file")]),
    (img2float,  realign,       [("out_file", "in_file")]),
    (extractref, realign,       [("roi_file", "ref_file")]),
    (realign,    plotrot,       [("par_file", "in_file")]),
    (realign,    plottrans,     [("par_file", "in_file")]),
    (realign,    plotdisp,      [("rms_files","in_file")]),
    (realign,    meanfunc1,     [("out_file", "in_file")]),
    (meanfunc1,  stripmean,     [("out_file", "in_file")]),
    (realign,    maskfunc1,     [("out_file", "in_file")]),
    (stripmean,  maskfunc1,     [("mask_file", "in_file2")]),
    (extractref, maskexample,   [("roi_file", "in_file")]),
    (stripmean,  maskexample,   [("mask_file", "in_file2")]),
    ])

# Register the fieldmap magnitude volume to the example func
regfieldmap = pe.MapNode(fsl.FLIRT(bins=256,
                                   cost="corratio",
                                   searchr_x=[-20, 20],
                                   searchr_y=[-20, 20],
                                   searchr_z=[-20, 20],
                                   dof=6,
                                   interp="trilinear"),
                          iterfield=["reference"],
                          name="regfieldmap")

# Resample the voxel-shift map to match the timeseries
resampvsm = pe.MapNode(interface=fsl.ApplyXfm(apply_xfm=True,
                                              interp="trilinear"),
                       iterfield=["reference", "in_matrix_file"],
                       name="resamplevsm")

# Dewarp the example func
dewarpexfunc = pe.MapNode(fsl.FUGUE(),
                          iterfield=["in_file", "shift_in_file", "mask_file"],
                          name="dewarpexfunc")

# Dewarp the timeseries
dewarpts = pe.MapNode(fsl.FUGUE(),
                          iterfield=["in_file", "shift_in_file", "mask_file"],
                          name="dewarptimeseries")


# Connect up the dewarping stages
preproc.connect([
    (inputnode,    regfieldmap,  [("fieldmap_brain", "in_file")]),
    (maskexample,  regfieldmap,  [("out_file", "reference")]),
    (regfieldmap,  resampvsm,    [("out_matrix_file", "in_matrix_file")]),
    (inputnode,    resampvsm,    [("voxel_shift_map", "in_file")]),
    (maskexample,  resampvsm,    [("out_file", "reference")]),
    (resampvsm,    dewarpexfunc, [("out_file", "shift_in_file")]),
    (maskexample,  dewarpexfunc, [("out_file", "in_file")]),
    (stripmean,    dewarpexfunc, [("mask_file", "mask_file")]),
    (resampvsm,    dewarpts,     [("out_file", "shift_in_file")]),
    (maskfunc1,    dewarpts,     [("out_file", "in_file")]),
    (stripmean,    dewarpts,     [("mask_file", "mask_file")]),
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

plotmean = pe.MapNode(fsl.PlotTimeSeries(title="Global Mean Intensity"),
                      iterfield=["in_file"],
                      name="plotmean")

# Make connections to ART
preproc.connect([
    (realign,    art,      [("par_file", "realignment_parameters")]),
    (maskfunc2,  art,      [("out_file", "realigned_files")]),
    (dilatemask, art,      [("out_file", "mask_file")]),
    (art,        plotmean, [("intensity_files", "in_file")]),
    ])

# Scale the median value each voxel in the run to 10000
meanscale = pe.MapNode(interface=fsl.ImageMaths(suffix="_gms"),
                     iterfield=["in_file","op_string"],
                     name="meanscale")

# High-pass filter the timeseries
highpass = pe.MapNode(interface=fsl.ImageMaths(suffix="_hpf",
                                               op_string="-bptf 64 -1"),
                      iterfield=["in_file"],
                      name="highpass")

# Get a new mean image from each functional run
meanfunc3 = pe.MapNode(interface=fsl.ImageMaths(op_string="-Tmean",
                                                suffix="_mean"),
                       iterfield=["in_file"],
                       name="meanfunc3")

# Slice the example func for reporting
meanslice = pe.MapNode(fsl.Slicer(image_width = 572,
                                  label_slices = False),
                       iterfield=["in_file"],
                       name="meanslice")
meanslice.inputs.sample_axial = 2

# Get the median value of the filtered timeseries
medianval2 = pe.MapNode(interface=fsl.ImageStats(op_string="-k %s -p 50"),
                        iterfield = ["in_file", "mask_file"],
                        name="medianval2")


# Functions to get fslmaths op strings
def getmeanscaleop(medianvals):
    """Get an fslmaths op string for intensity normalization."""
    return ["-mul %.10f"%(10000./val) for val in medianvals]

def gethpfop(cutoff):
    """Get an fslmaths op string for high-pass filtering given a cutoff in TRs."""
    return "-bptf %d -1"%(cutoff/2)

# Make connections to the intensity normalization and temporal filtering
preproc.connect([
    (maskfunc2, meanscale,  [("out_file", "in_file")]),
    (medianval, meanscale,  [(("out_stat", getmeanscaleop), "op_string")]),
    (inputnode, highpass,   [(("hpf_cutoff", gethpfop), "op_string")]),
    (meanscale, highpass,   [("out_file", "in_file")]),
    (highpass,  meanfunc3,  [("out_file", "in_file")]),
    (meanfunc3, meanslice,  [("out_file", "in_file")]),
    (highpass,  medianval2, [("out_file", "in_file")]),
    (threshold, medianval2, [("out_file", "mask_file")]),
    ])

# Merge the median values with the mean functional images into a coupled list
mergenode = pe.Node(interface=util.Merge(2, axis="hstack"),
                    name="merge")

# Smooth in the volume with SUSAN
smooth = pe.MapNode(interface=fsl.SUSAN(),
                    iterfield=["in_file", "brightness_threshold", "usans"],
                    name="smooth")

# Mask the smoothed data with the dilated mask
masksmoothfunc = pe.MapNode(interface=fsl.ImageMaths(suffix="_mask",
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
    (meanfunc3,  mergenode,      [("out_file", "in1")]),
    (medianval2, mergenode,      [("out_stat", "in2")]),
    (highpass,   smooth,         [("out_file", "in_file")]),
    (inputnode,  smooth,         [("smooth_fwhm", "fwhm")]),
    (medianval2, smooth,         [(("out_stat", getbtthresh), "brightness_threshold")]),
    (mergenode,  smooth,         [(("out", getusans), "usans")]),
    (smooth,     masksmoothfunc, [("smoothed_file", "in_file")]),
    (dilatemask, masksmoothfunc, [("out_file", "in_file2")]),
    ])

# Define the outputs of the preprocessing that will be used by the model
outputnode = pe.Node(util.IdentityInterface(fields=["smoothed_timeseries",
                                                    "unsmoothed_timeseries",
                                                    "example_func",
                                                    "mean_func",
                                                    "functional_mask",
                                                    "realignment_parameters",
                                                    "outlier_files"]),
                     name="outputspec")

# Define the outputs that will go into a report
reportnode = pe.Node(util.IdentityInterface(fields=["example_func",
                                                    "mean_func",
                                                    "intensity_plot",
                                                    "outlier_volumes",
                                                    "rotation_plot",
                                                    "translation_plot",
                                                    "displacement_plot"]),
                     name="report")

# Make connections to the output nodes
preproc.connect([
    (highpass,       outputnode, [("out_file", "unsmoothed_timeseries")]),
    (masksmoothfunc, outputnode, [("out_file", "smoothed_timeseries")]),
    (extractref,     outputnode, [("roi_file", "example_func")]),
    (meanfunc2,      outputnode, [("out_file", "mean_func")]),
    (dilatemask,     outputnode, [("out_file", "functional_mask")]),
    (realign,        outputnode, [("par_file", "realignment_parameters")]),
    (art,            outputnode, [("outlier_files", "outlier_files")]),
    (art,            reportnode, [("outlier_files", "outlier_volumes")]),
    (plotmean,       reportnode, [("out_file", "intensity_plot")]),
    (exampleslice,   reportnode, [("out_file", "example_func")]),
    (meanslice,      reportnode, [("out_file", "mean_func")]),
    (plotrot,        reportnode, [("out_file", "rotation_plot")]),
    (plottrans,      reportnode, [("out_file", "translation_plot")]),
    (plotdisp,       reportnode, [("out_file", "displacement_plot")]),
    ])

