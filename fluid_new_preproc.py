"""
    Preprocessing module for Fluid Intelligence fMRI paradigms.
"""


import os                               
import sys
import nipype.interfaces.io as nio       
import nipype.interfaces.fsl as fsl       
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util   
import nipype.pipeline.engine as pe         
import nipype.algorithms.rapidart as ra      
import nibabel


preproc = pe.Workflow(name="preproc")

"""
Set up a node to define all inputs required for the preprocessing workflow
"""

inputnode = pe.Node(interface=util.IdentityInterface(fields=["func",
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

extract_ref = pe.MapNode(interface=fsl.ExtractROI(t_size=1,
                                                  t_min=0),
                         iterfield=["in_file"],
                         name = "extractref")

preproc.connect(img2float, "out_file", extract_ref, "in_file")

"""
Realign each functional run to the middle volume of the run
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
Generate linear registration matricies for each run to the structural image
"""

reg2struct = pe.MapNode(interface=fs.BBRegister(contrast_type="t2",
                                                init = "fsl"),
                        iterfield=["source_file"],
                        name="reg2struct")

preproc.connect(meanfunc, "out_file", reg2struct, "source_file")
preproc.connect(inputnode, "subject_id", reg2struct, "subject_id")

"""
Apply this transformation to each mean and full functional timeseires
"""

meanxfm = pe.MapNode(interface=fs.ApplyVolTransform(fs_target=True,
                                                    no_resample=True),
                     iterfield=["reg_file","source_file"],
                     name="applyxfm2mean")

preproc.connect(meanfunc, "out_file", meanxfm, "source_file")
preproc.connect(reg2struct, "out_reg_file", meanxfm, "reg_file")

tsxfm = pe.MapNode(interface=fs.ApplyVolTransform(fs_target=True,
                                                  no_resample=True),
                   iterfield=["reg_file","source_file"],
                   name="applyxfm2ts")

preproc.connect(motion_correct, "out_file", tsxfm, "source_file")
preproc.connect(reg2struct, "out_reg_file", tsxfm, "reg_file")

"""
Strip the skull from the transformed mean functional from the first run
"""

def pickfirst(files):
    """Return the first file from a list"""
    if isinstance(files, list):
        return files[0]
    else:
        return files

meanfuncmask = pe.Node(interface=fsl.BET(mask = True,
                                         no_output=True,
                                         frac = 0.3),
                       name = "meanfuncmask")
preproc.connect(meanxfm, ("transformed_file", pickfirst), meanfuncmask, "in_file")

"""
Mask the functional runs with the extracted mask
"""

maskfunc = pe.MapNode(interface=fsl.ImageMaths(suffix="_bet",
                                               op_string="-mas"),
                      iterfield=["in_file"],
                      name = "maskfunc")
preproc.connect(tsxfm, "transformed_file", maskfunc, "in_file")
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

threshold = pe.Node(interface=fsl.ImageMaths(out_data_type="char",
                                             suffix="_thresh"),
                       name="threshold")
preproc.connect(maskfunc, ("out_file", pickfirst), threshold, "in_file")

"""
Define a function to get 10% of the intensity
"""

def getthreshop(thresh):
    return "-thr %.10f -Tmin -bin"%(0.1*thresh[0][1])
preproc.connect(getthresh, ("out_stat", getthreshop), threshold, "op_string")

"""
Determine the median value of the functional runs using the mask
"""

medianval = pe.MapNode(interface=fsl.ImageStats(op_string="-k %s -p 50"),
                       iterfield = ["in_file"],
                       name="medianval")
preproc.connect(motion_correct, "out_file", medianval, "in_file")
preproc.connect(threshold, "out_file", medianval, "mask_file")

"""
Dilate the mask
"""

dilatemask = pe.Node(interface=fsl.ImageMaths(suffix="_dil",
                                              op_string="-dilF"),
                       name="dilatemask")
preproc.connect(threshold, "out_file", dilatemask, "in_file")

"""
Mask the motion corrected functional runs with the dilated mask
"""

maskfunc2 = pe.MapNode(interface=fsl.ImageMaths(suffix="_mask",
                                                op_string="-mas"),
                      iterfield=["in_file"],
                      name="maskfunc2")
preproc.connect(motion_correct, "out_file", maskfunc2, "in_file")
preproc.connect(dilatemask, "out_file", maskfunc2, "in_file2")

"""
Determine the mean image from each functional run
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
Identitynode to set fwhm
"""
smoothval = pe.Node(interface=util.IdentityInterface(fields=["fwhm"]),
                    name="smoothval")
                    
                       
"""
Smooth each run using SUSAN with the brightness threshold set to 75% of the
median value for each run and a mask consituting the mean functional
"""

smooth = pe.MapNode(interface=fsl.SUSAN(),
                    iterfield=["in_file", "brightness_threshold","usans"],
                    name="smooth")

"""
Define a function to get the brightness threshold for SUSAN
"""

def getbtthresh(medianvals):
    return [0.75*val for val in medianvals]
preproc.connect(smoothval, "fwhm", smooth, "fwhm")
preproc.connect(maskfunc2, "out_file", smooth, "in_file")
preproc.connect(medianval, ("out_stat", getbtthresh), smooth, "brightness_threshold")
preproc.connect(mergenode, 
    ("out", lambda x: [[tuple([val[0],0.75*val[1]])] for val in x]), smooth, "usans")

"""
Mask the smoothed data with the dilated mask
"""

maskfunc3 = pe.MapNode(interface=fsl.ImageMaths(suffix="_mask",
                                                op_string="-mas"),
                      iterfield=["in_file"],
                      name="maskfunc3")
preproc.connect(smooth, "smoothed_file", maskfunc3, "in_file")
preproc.connect(dilatemask, "out_file", maskfunc3, "in_file2")


concatnode = pe.Node(interface=util.Merge(2),
                     name="concat")
preproc.connect(maskfunc2,("out_file", lambda x:[x]), concatnode, "in1")
preproc.connect(maskfunc3,("out_file", lambda x:[x]), concatnode, "in2")

selectnode = pe.Node(interface=util.Select(),name="select")

preproc.connect(concatnode, "out", selectnode, "inlist")

def chooseindex(fwhm):
    if fwhm<1:
        return [0]
    else:
        return [1]
    
preproc.connect(smoothval, ("fwhm", chooseindex), selectnode, "index")


"""
Scale the median value of the run is set to 10000
"""

meanscale = pe.MapNode(interface=fsl.ImageMaths(suffix="_gms"),
                      iterfield=["in_file","op_string"],
                      name="meanscale")
preproc.connect(selectnode, "out", meanscale, "in_file")

def getmeanscale(medianvals):
    """Get the scaling factor for intensity normalization"""
    return ["-mul %.10f"%(10000./val) for val in medianvals]
preproc.connect(medianval, ("out_stat", getmeanscale), meanscale, "op_string")

"""
Perform temporal highpass filtering on the data
"""

highpass = pe.MapNode(interface=fsl.ImageMaths(suffix="_tempfilt"),
                      iterfield=["in_file"],
                      name="highpass")
preproc.connect(meanscale, "out_file", highpass, "in_file")

"""
Generate a mean functional image from the first run
"""

meanfunc3 = pe.MapNode(interface=fsl.ImageMaths(op_string="-Tmean",
                                                suffix="_mean"),
                       iterfield=["in_file"],
                      name="meanfunc3")
preproc.connect(highpass, ("out_file", pickfirst), meanfunc3, "in_file")

"""
Use :class:`nipype.algorithms.rapidart` to determine which of the
images in the functional series are outliers based on deviations in
intensity and/or movement.
"""

art = pe.Node(interface=ra.ArtifactDetect(use_differences = [True, False],
                                          use_norm = True,
                                          zintensity_threshold = 3,
                                          norm_threshold = 1,
                                          parameter_source = "FSL",
                                          mask_type = "file"),
              name="art")


preproc.connect([(motion_correct, art, [("par_file","realignment_parameters")]),
                 (maskfunc2, art, [("out_file","realigned_files")]),
                 (dilatemask, art, [("out_file", "mask_file")]),
                 ])


