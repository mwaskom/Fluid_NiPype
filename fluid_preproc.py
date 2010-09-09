# vi: set ft=python sts=4 ts=4 sw=4 et:
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
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

inputnode = pe.Node(interface=util.IdentityInterface(fields=['func',
                                                             'struct']),
                    name='inputspec')

"""
Convert functional images to float representation. Since there can be more than
one functional run we use a MapNode to convert each run.
"""

img2float = pe.MapNode(interface=fsl.ImageMaths(out_data_type='float',
                                             op_string = '',
                                             suffix='_dtype'),
                       iterfield=['in_file'],
                       name='img2float')
preproc.connect(inputnode, 'func', img2float, 'in_file')

"""
Extract the middle volume of the first run as the reference
"""

extract_ref = pe.MapNode(interface=fsl.ExtractROI(t_size=1),
                         iterfield=["in_file"],
                         name = 'extractref')

"""
Define a function to pick the first file from a list of files
"""

def pickfirst(files):
    if isinstance(files, list):
        return files[0]
    else:
        return files

preproc.connect(img2float, 'out_file', extract_ref, 'in_file')

"""
Define a function to return the 1 based index of the middle volume
"""

def getmiddlevolume(func):
    funcfile = func
    if isinstance(func, list):
        funcfile = func[0]
    _,_,_,timepoints = nibabel.load(funcfile).get_shape()
    return (timepoints/2)-1

preproc.connect(inputnode, ('func', getmiddlevolume), extract_ref, 't_min')

"""
Realign the functional runs to the middle volume of that run
"""

motion_correct = pe.MapNode(interface=fsl.MCFLIRT(save_mats = True,
                                                  save_plots = True,
                                                  save_rms = True),
                            name='realign',
                            iterfield = ['in_file', "ref_file"])
preproc.connect(img2float, 'out_file', motion_correct, 'in_file')
preproc.connect(extract_ref, 'roi_file', motion_correct, 'ref_file')

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
Get a mean image for each realigned run
"""

meanfunc = pe.MapNode(interface=fsl.ImageMaths(op_string='-Tmean',
                                               suffix='_mean'),
                       iterfield=['in_file'],
                       name='meanfunc')

preproc.connect(realign, "out_file", meanfunc, "in_file")

"""
Generate linear registration matricies for each run to the anatomical
volume and apply the registration to the timeseries(es)
"""

reg2struct = pe.MapNode(interface=fs.BBRegister(contrast_type="bold",
                                            init = "fsl"),
                        iterfield=["source_file"],
                        name="reg2struct")

preproc.connect(meanfunc, "out_file", reg2struct, "source_file")

"""
Get the brainmask volume from recon-all
"""

masksource = pe.Node(interface=nio.FreeSurferSource(),
                     name="masksource")

"""
Resample the brainmask image into native functional space
"""

resamplemask = pe.Node(interface=fs.ApplyVolTransform(inverse=True),
                       name="resamplemask")

preproc.connect([(masksource, resamplemask, [("brainmask", "source_file")]),
                  (reg2struct, resamplemask, [("out_reg_file", "reg_file")])])

"""
Binarize the brainmask image in functional space to create a true mask
"""

binmask = pe.Node(interface=fs.Binarize(min=10),
                  name="binarizemask")

preproc.connect(resamplemask, "out_file", binmask, "in_file")

"""
Convert the mask image from mgz to nifti so it can be used by FSL tools
"""

niftimask = pe.Node(interface=fs.MRIConvert(out_type="niigz"),
                    name="niftimask")

preproc.connect(binmask, "out_file", niftimask, "in_file")


"""
Use this mask image to skullstrip each realigned functional run
"""

skullstrip = pe.MapNode(interface=fsl.ImageMaths(suffix="_strip",
                                                 op_string="-mas"),
                        iterfield=["in_file"],
                        name="skullstrip")

preproc.connect([(motion_correct, skullstrip, [('out_file', 'in_file')]),
                 (niftimask, skullstrip, [("out_file", "in_file2")])
                 ]) 

"""
Determine the 2nd and 98th percentile intensities of each functional run
"""

getthresh = pe.MapNode(interface=fsl.ImageStats(op_string='-p 2 -p 98'),
                       iterfield = ['in_file'],
                       name='getthreshold')
preproc.connect(skullstrip, 'out_file', getthresh, 'in_file')

"""
Threshold the first run of the functional data at 10% of the 98th percentile
"""

threshold = pe.Node(interface=fsl.ImageMaths(out_data_type='char',
                                             suffix='_thresh'),
                       name='threshold')
preproc.connect(skullstrip, ('out_file', pickfirst), threshold, 'in_file')


def getthreshop(thresh):
    """Return 10% of the intensity as an fslmaths option string"""
    return '-thr %.10f -Tmin -bin'%(0.1*thresh[0][1])
preproc.connect(getthresh, ('out_stat', getthreshop), threshold, 'op_string')

"""
Determine the median value of the functional runs using the mask
"""

medianval = pe.MapNode(interface=fsl.ImageStats(op_string='-k %s -p 50'),
                       iterfield = ['in_file'],
                       name='medianval')
preproc.connect(motion_correct, 'out_file', medianval, 'in_file')
preproc.connect(threshold, 'out_file', medianval, 'mask_file')

"""
Merge the median values with the mean functional images into a coupled list
"""

mergenode = pe.Node(interface=util.Merge(2, axis='hstack'),
                    name='merge')
preproc.connect(coregister, 'out_file', mergenode, 'in1')
preproc.connect(medianval,'out_stat', mergenode, 'in2')

                       
"""
Smooth each run using SUSAN with the brightness threshold set to 75% of the
median value for each run and a mask consituting the mean functional
"""

smooth = pe.MapNode(interface=fsl.SUSAN(),
                    iterfield=['in_file', 'brightness_threshold','usans'],
                    name='smooth')

"""
Define a function to get the brightness threshold for SUSAN
"""

def getbtthresh(medianvals):
    return [0.75*val for val in medianvals]

preproc.connect(applyxfm, 'out_file', smooth, 'in_file')
preproc.connect(medianval, ('out_stat', getbtthresh), smooth, 'brightness_threshold')
preproc.connect(mergenode, ('out', lambda x: [[tuple([val[0],0.75*val[1]])] for val in x]), smooth, 'usans')

"""
Mask the smoothed data with the dilated mask
"""

maskfunc3 = pe.MapNode(interface=fsl.ImageMaths(suffix='_mask',
                                                op_string='-mas'),
                      iterfield=['in_file'],
                      name='maskfunc3')
preproc.connect(smooth, 'smoothed_file', maskfunc3, 'in_file')
preproc.connect(dilatemask, 'out_file', maskfunc3, 'in_file2')

"""
Scale each volume of the run so that the median value of the run is set to 10000
"""

intnorm = pe.MapNode(interface=fsl.ImageMaths(suffix='_intnorm'),
                      iterfield=['in_file','op_string'],
                      name='intnorm')
preproc.connect(maskfunc3, 'out_file', intnorm, 'in_file')

"""
Define a function to get the scaling factor for intensity normalization
"""

def getinormscale(medianvals):
    return ['-mul %.10f'%(10000./val) for val in medianvals]
preproc.connect(medianval, ('out_stat', getinormscale), intnorm, 'op_string')

"""
Perform temporal highpass filtering on the data
"""

highpass = pe.MapNode(interface=fsl.ImageMaths(suffix='_tempfilt'),
                      iterfield=['in_file'],
                      name='highpass')
preproc.connect(intnorm, 'out_file', highpass, 'in_file')

"""
Convert the intensity normalization output to .nii so SPM can read it
"""

convert = pe.MapNode(interface=fs.MRIConvert(out_type="nii"), 
                     iterfield = "in_file",
                     name="unzip_intnorm")
preproc.connect(intnorm, "out_file", convert, "in_file")

"""
Generate a mean functional image from the first run
"""

meanfunc3 = pe.MapNode(interface=fsl.ImageMaths(op_string='-Tmean',
                                                suffix='_mean'),
                       iterfield=['in_file'],
                      name='meanfunc3')
preproc.connect(highpass, ('out_file', pickfirst), meanfunc3, 'in_file')


"""
Use rapidart to determine which of the images in the functional series are 
outliers based on deviations in intensity and/or movement.
"""

art = pe.Node(interface=ra.ArtifactDetect(use_differences = [False,True],
                                          use_norm = True,
                                          norm_threshold = 0.5,
                                          zintensity_threshold = 3,
                                          parameter_source = 'FSL',
                                          mask_type = 'file'),
              name="art")


preproc.connect([(motion_correct, art, [('par_file','realignment_parameters')]),
                 (skullstrip, art, [('out_file','realigned_files')]),
                 (niftimask, art, [('out_file', 'mask_file')]),
                 ])
