"""
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

"""
Set up a node to define all inputs required for the preprocessing workflow
"""

inputnode = pe.Node(interface=util.IdentityInterface(fields=["func",
                                                             "ep2d_t1w",
                                                             "subject_id"]),
                    name="inputspec")

"""
Convert functional images to float representation. Since there can be more than
one functional run we use a MapNode to convert each run.
"""

img2float = pe.MapNode(interface=fsl.ImageMaths(out_data_type="float",
                                             op_string = "",
                                             suffix="_dtype"),
                       iterfield=["in_file"],
                       name="img2float")
preproc.connect(inputnode, "func", img2float, "in_file")

"""
Extract the middle volume of the first run as the reference
"""


def getmiddlevolume(func):
    """Return the middle volume index to use as a reference in MCFLIRT"""
    funcfile = func
    if isinstance(func, list):
        funcfile = func[0]
    _,_,_,timepoints = nib.load(funcfile).get_shape()
    return (timepoints/2)-1

extract_ref = pe.MapNode(interface=fsl.ExtractROI(t_size=1),
                         iterfield=["in_file"],
                         name = "extractref")

preproc.connect(img2float, "out_file", extract_ref, "in_file")
preproc.connect(inputnode, ("func", getmiddlevolume), extract_ref, "t_min")

"""
Realign each functional run to the middle volume of that run
"""

motion_correct = pe.MapNode(interface=fsl.MCFLIRT(save_mats = True,
                                                  save_plots = True,
                                                  save_rms = True),
                            name="realign",
                            iterfield = ["in_file", "ref_file"])
preproc.connect(img2float, "out_file", motion_correct, "in_file")
preproc.connect(extract_ref, "roi_file", motion_correct, "ref_file")

"""
Plot the motion and displacement estimates from MCFLIRT
"""

plotmotion = pe.MapNode(interface=fsl.PlotMotionParams(in_source="fsl"),
                        name="plotmotion",
                        iterfield=["in_file"])
plotmotion.iterables = ("plot_type", ["rotations","translations"])        

preproc.connect(motion_correct, "par_file", plotmotion, "in_file")


plotdisp = pe.MapNode(interface=fsl.PlotMotionParams(in_source="fsl",
                                                     plot_type="displacement"),
                      name="plotdisplacement",
                      iterfield=["in_file"])

preproc.connect(motion_correct, "rms_files", plotdisp, "in_file")

"""
Extract the mean volume of each functional run
"""

meanfunc = pe.MapNode(interface=fsl.ImageMaths(op_string = "-Tmean",
                                               suffix="_mean"),
                      iterfield=["in_file"],
                      name="meanfunc")
preproc.connect(motion_correct, "out_file", meanfunc, "in_file")

"""
Generate linear registration matrices for each run to the t1-weighted EPI
"""

reg2struct = pe.MapNode(fsl.FLIRT(),
                        iterfield=["in_file"],
                        name="reg2struct")

preproc.connect(inputnode, "ep2d_t1w", reg2struct, "reference")
preproc.connect(meanfunc, "out_file", reg2struct, "in_file")



"""
DEPRECATED

reg2struct = pe.MapNode(interface=fs.BBRegister(contrast_type="t2",
                                                init = "fsl"),
                        iterfield=["source_file"],
                        name="reg2struct")

preproc.connect(meanfunc, "out_file", reg2struct, "source_file")
preproc.connect(inputnode, "subject_id", reg2struct, "subject_id")

meanxfm = pe.MapNode(interface=fs.ApplyVolTransform(fs_target=True,
                                                    no_resample=True),
                     iterfield=["reg_file","source_file"],
                     name="applyxfm2mean")

preproc.connect(meanfunc, "out_file", meanxfm, "source_file")
preproc.connect(reg2struct, "out_reg_file", meanxfm, "reg_file")
meanmean = pe.Node(interface=fs.Concatenate(stats="mean"),
                   name="meanmean")
                   
preproc.connect(meanxfm, "transformed_file", meanmean, "in_files")


surfreg = pe.Node(interface=fs.BBRegister(contrast_type="t2",
                                          init = "fsl"),
                  name="surfregister")

preproc.connect(meanmean, "concatenated_file", surfreg, "source_file")
preproc.connect(inputnode, "subject_id", surfreg, "subject_id")

"""

"""
Register the t1 weighted target image to the original anatomical
"""

surfreg = pe.Node(fs.BBRegister(init="fsl", contrast_type="t1"),
                  name="surfreg")

preproc.connect(inputnode, "subject_id", surfreg, "subject_id")
preproc.connect(inputnode, "ep2d_t1w", surfreg, "source_file")

"""
Get the brainmask produced by recon-all
"""

fssource = pe.Node(interface=nio.FreeSurferSource(), name="fssource")

preproc.connect(inputnode, "subject_id", fssource, "subject_id")

"""
Invert the surfreg transform to put the brainmask into functional space
"""

maskxfm = pe.Node(interface=fs.ApplyVolTransform(inverse=True),
                  name="maskxfm")

preproc.connect([(fssource, maskxfm, [("brainmask", "target_file")]),
                 (inputnode, maskxfm, [("ep2d_t1w", "source_file")]),
                 (surfreg, maskxfm, [("out_reg_file", "reg_file")])
                 ])

"""
Binarize the brainmask in functional space
"""

threshmask = pe.Node(interface=fs.Binarize(min=10), name = "threshmask")

preproc.connect(maskxfm, "transformed_file", threshmask, "in_file")

"""
Convert the mask to nifti so FSL tools can read it
"""

mask2nii = pe.Node(interface=fs.MRIConvert(out_type="niigz"), name="mask2nii")

preproc.connect(threshmask, "binary_file", mask2nii, "in_file")

"""
Strip the skull from the mean functional from each run
"""

meanfuncmask = pe.MapNode(interface=fsl.BET(mask = True,
                                         no_output=True,
                                         frac = 0.3),
                          iterfield = "in_file",
                          name = "meanfuncmask")
preproc.connect(meanfunc, "out_file", meanfuncmask, "in_file")

"""
Mask the functional runs with the extracted mask
"""

maskfunc = pe.MapNode(interface=fsl.ImageMaths(suffix="_bet",
                                               op_string="-mas"),
                      iterfield=["in_file", "in_file2"],
                      name = "maskfunc")
preproc.connect(motion_correct, "out_file", maskfunc, "in_file")
preproc.connect(meanfuncmask, "mask_file", maskfunc, "in_file2")

"""
Determine the 2nd and 98th percentile intensities of each functional run
"""

getthresh = pe.MapNode(interface=fsl.ImageStats(op_string="-p 2 -p 98"),
                       iterfield = ["in_file"],
                       name="getthreshold")
preproc.connect(maskfunc, "out_file", getthresh, "in_file")

"""
Threshold the first run of the functional data at 10% of the 98th percentile
"""

threshold = pe.MapNode(interface=fsl.ImageMaths(out_data_type="char",
                                                suffix="_thresh"),
                       iterfield = ["in_file"],
                       name="threshold")
preproc.connect(maskfunc, "out_file", threshold, "in_file")

def getthreshop(thresh):
    """Return 10% of the intensity"""
    return "-thr %.10f -Tmin -bin"%(0.1*thresh[0][1])
preproc.connect(getthresh, ("out_stat", getthreshop), threshold, "op_string")

"""
Determine the median value of the functional runs using the mask
"""

medianval = pe.MapNode(interface=fsl.ImageStats(op_string="-k %s -p 50"),
                       iterfield = ["in_file", "mask_file"],
                       name="medianval")
preproc.connect(motion_correct, "out_file", medianval, "in_file")
preproc.connect(threshold, "out_file", medianval, "mask_file")

"""
Dilate the mask
"""

dilatemask = pe.MapNode(interface=fsl.ImageMaths(suffix="_dil",
                                                 op_string="-dilF"),
                        iterfield=["in_file"],
                        name="dilatemask")
preproc.connect(threshold, "out_file", dilatemask, "in_file")

"""
Mask the motion corrected functional runs with the dilated mask
"""

maskfunc2 = pe.MapNode(interface=fsl.ImageMaths(suffix="_mask",
                                                op_string="-mas"),
                      iterfield=["in_file", "in_file2"],
                      name="maskfunc2")
preproc.connect(motion_correct, "out_file", maskfunc2, "in_file")
preproc.connect(dilatemask, "out_file", maskfunc2, "in_file2")

"""
Use RapidART to determine which of the
images in the functional series are outliers based on deviations in
intensity and/or movement.
"""

art = pe.MapNode(interface=ra.ArtifactDetect(use_differences = [True, False],
                                             use_norm = True,
                                             zintensity_threshold = 3,
                                             norm_threshold = 1,
                                             parameter_source = "FSL",
                                             mask_type = "file"),
                 iterfield=["realignment_parameters","realigned_files","mask_file"],
                 name="art")


preproc.connect([(motion_correct, art, [("par_file","realignment_parameters")]),
                 (maskfunc2, art, [("out_file","realigned_files")]),
                 (dilatemask, art, [("out_file", "mask_file")]),
                 ])
                 
"""
Determine a new mean image from each functional run
"""

meanfunc2 = pe.MapNode(interface=fsl.ImageMaths(op_string="-Tmean",
                                                suffix="_mean"),
                       iterfield=["in_file"],
                       name="meanfunc2")
preproc.connect(maskfunc2, "out_file", meanfunc2, "in_file")

"""
Merge the median values with the mean functional images into a coupled list
"""

mergenode = pe.Node(interface=util.Merge(2, axis="hstack"),
                    name="merge")
preproc.connect(meanfunc2,"out_file", mergenode, "in1")
preproc.connect(medianval,"out_stat", mergenode, "in2")


"""
Identitynodes to set fwhm for smoothing
"""
volsmoothval = pe.Node(interface=util.IdentityInterface(fields=["fwhm"]),
                       name="volsmoothval")
surfsmoothval = pe.Node(interface=util.IdentityInterface(fields=["fwhm"]),
                        name="surfsmoothval")
                       
"""
Smooth each run in the volume using SUSAN with the brightness threshold set to
75% of the median value for each run and a mask consituting the mean functional
"""

volsmooth = pe.MapNode(interface=fsl.SUSAN(),
                       iterfield=["in_file", "brightness_threshold", "usans"],
                       name="volsmooth")


def getbtthresh(medianvals):
    """Get the brightness threshold for SUSAN"""
    return [0.75*val for val in medianvals]
preproc.connect(volsmoothval, "fwhm", volsmooth, "fwhm")
preproc.connect(maskfunc2, "out_file", volsmooth, "in_file")
preproc.connect(medianval, ("out_stat", getbtthresh), volsmooth, "brightness_threshold")
preproc.connect(mergenode, 
    ("out", lambda x: [[tuple([val[0],0.75*val[1]])] for val in x]), volsmooth, "usans")

"""
Use Freesurfer's vol/surf smoothing to smooth the images for surface analysis
"""

surfsmooth = pe.MapNode(interface=fs.Smooth(proj_frac_avg=(0,1,0.1)),
                        iterfield=["in_file","reg_file"],
                        name="surfsmooth")

preproc.connect(maskfunc2, "out_file", surfsmooth, "in_file")
preproc.connect(surfsmoothval, "fwhm", surfsmooth, "surface_fwhm")
preproc.connect(reg2struct, "out_file", surfsmooth, "reg_file")

"""
Mask the smoothed data with the dilated mask
"""

volmaskfunc3 = pe.MapNode(interface=fsl.ImageMaths(suffix="_mask",
                                                op_string="-mas"),
                      iterfield=["in_file","in_file2"],
                      name="volmaskfunc3")
preproc.connect(volsmooth, "smoothed_file", volmaskfunc3, "in_file")
preproc.connect(dilatemask, "out_file", volmaskfunc3, "in_file2")

"""
Scale the median value of the run is set to 10000
"""

volmeanscale = pe.MapNode(interface=fsl.ImageMaths(suffix="_gms"),
                          iterfield=["in_file","op_string"],
                          name="volmeanscale")
preproc.connect(volsmooth, "smoothed_file", volmeanscale, "in_file")

def getmeanscale(medianvals):
    """Get the scaling factor for intensity normalization"""
    return ["-mul %.10f"%(10000./val) for val in medianvals]
preproc.connect(medianval, ("out_stat", getmeanscale), volmeanscale, "op_string")

surfmeanscale = pe.MapNode(interface=fsl.ImageMaths(suffix="_gms"),
                           iterfield=["in_file","op_string"],
                           name="surfmeanscale")
preproc.connect(surfsmooth, "smoothed_file", surfmeanscale, "in_file")
preproc.connect(medianval, ("out_stat", getmeanscale), surfmeanscale, "op_string")

"""
Perform temporal highpass filtering on the data
"""

volhighpass = pe.MapNode(interface=fsl.ImageMaths(suffix="_vol_hpf"),
                         iterfield=["in_file"],
                         name="volhighpass")
preproc.connect(volmeanscale, "out_file", volhighpass, "in_file")

surfhighpass = pe.MapNode(interface=fsl.ImageMaths(suffix="_surf_hpf"),
                          iterfield=["in_file"],
                          name="surfhighpass")
preproc.connect(surfmeanscale, "out_file", surfhighpass, "in_file")

"""
Use ART to find stimulus correlated motion
"""

unzip = pe.MapNode(interface=fs.MRIConvert(out_type="nii"),
                   iterfield=["in_file"],
                   name="unzipforart")

preproc.connect(volhighpass, "out_file", unzip, "in_file")

artspec = pe.Node(interface=model.SpecifyModel(concatenate_runs=False),
                  name="artspec")

preproc.connect(unzip, "out_file", artspec, "functional_runs")

artdesign = pe.Node(interface=spm.Level1Design(), 
                    name="artdesign")

preproc.connect(artspec, "session_info", artdesign, "session_info")

stimcorr = pe.Node(interface=ra.StimulusCorrelation(concatenated_design=False),
                   name="stimcorr")

preproc.connect(artdesign, "spm_mat_file", stimcorr, "spm_mat_file")

preproc.connect([
    (motion_correct, artspec, [("par_file", "realignment_parameters")]),
    (art, artspec, [("outlier_files", "outlier_files")]),
    (art, stimcorr, [("intensity_files", "intensity_values")]),
    ])

