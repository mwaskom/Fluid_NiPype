# vi: set ft=python sts=4 ts=4 sw=4 et:
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
"""
    Nipype interface script for the Fluid Intelligence project.
"""

import os
import re
import sys
import shutil
from copy import deepcopy
from datetime import datetime

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.spm as spm
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util

# Catch when we"re not using the right environment
if (not pe.__file__.startswith("/software/python/nipype0.3") 
    and not pe.__file__.startswith("/u2/mwaskom/nipype")):
    sys.exit("ERROR: Not using nipype0.3")

import argparse

from fluid_preproc import preproc
from fluid_fsl_model import fsl_vol_model, fsl_surf_model
from fluid_spm_model import spm_modelfit
from fluid_fixed_fx import vol_fixed_fx, surf_fixed_fx
from fluid_source import infosource, datasource

import fluid_utility_funcs as fuf

""" Handle command line arguments to control the analysis """
parser = argparse.ArgumentParser(description="Main interface for GFluid NiPype code.")
parser.add_argument("-paradigm", metavar="paradigm", required=True,
                    help="experimental paradigm")
parser.add_argument("-subject", dest="subjects", 
                    metavar="subject_id", action="append",
                    help="run pypeline for subject(s)")
parser.add_argument("-norun",dest="run",action="store_false",
                    help="do not run the pypeline")
parser.add_argument("-nograph",dest="write_graph",action="store_false",
                    help="do not write a graph")
parser.add_argument("-inseries",action="store_true",
                    help="force running in series")
parser.add_argument("-runspm",action="store_true",
                    help="run spm first level")
args = parser.parse_args()
""" Dynamically import the experiment file """
exp = __import__("%s_experiment" % args.paradigm)

""" Subjects.  This won"t stay hardcorded like this """
if hasattr(exp, "subject_list"):
    subject_list = exp.subject_list
else:
    subject_list = ["gf%02d"%id for id in [5, 9, 13, 14]]
if hasattr(exp, "exclude_subjects"):    
    subject_list = [subj for subj in subject_list if subj not in exp.exclude_subjects]

if args.subjects:
    subject_list = args.subjects

""" Set experimental source attributes. """
datasource.inputs.template = exp.source_template
datasource.inputs.template_args = exp.template_args

exp.data_dir = datasource.inputs.base_directory

""" Define the level 1 pipeline"""
firstlevel = pe.Workflow(name= "level1")

""" Tell infosource to iterate over all subjects """
infosource.iterables = ("subject_id", subject_list)

""" Add experiment-sepcific information to the workflows """
highpass_opstring = "-bptf %s -1" % (exp.hpcutoff/exp.TR)
preproc.inputs.volhighpass.op_string = highpass_opstring
preproc.inputs.surfhighpass.op_string = highpass_opstring

surfsmoothvalnode = preproc.get_node("surfsmoothval")
#surfsmoothvalnode.iterables = ("surface_fwhm", [0.,5.])
preproc.inputs.surfsmoothval.fwhm = 5.
preproc.inputs.volsmoothval.fwhm = 5.

""" Connect the workflows """
firstlevel.connect([(infosource, datasource, 
                        [("subject_id", "subject_id")]),
                    (infosource, preproc,
                        [("subject_id", "inputspec.subject_id")]),
                    (datasource, preproc, 
                        [("func", "inputspec.func")])
                   ])

contrasts = exp.contrasts

fsl_vol_model.inputs.modelspec.time_repetition = exp.TR
fsl_vol_model.inputs.modelspec.high_pass_filter_cutoff = exp.hpcutoff
fsl_vol_model.inputs.modelspec.input_units = exp.units
fsl_vol_model.inputs.modelspec.output_units = exp.units

fsl_vol_model.inputs.level1design.contrasts = contrasts
fsl_vol_model.inputs.level1design.interscan_interval = exp.TR
fsl_vol_model.inputs.level1design.bases = exp.fsl_bases

selectvolcontrast = fsl_vol_model.get_node("selectcontrast")
selectvolcontrast.iterables = ('index', [[i] for i in range(len(contrasts))])

fsl_surf_model.inputs.modelspec.time_repetition = exp.TR
fsl_surf_model.inputs.modelspec.high_pass_filter_cutoff = exp.hpcutoff
fsl_surf_model.inputs.modelspec.input_units = exp.units
fsl_surf_model.inputs.modelspec.output_units = exp.units

fsl_surf_model.inputs.level1design.contrasts = contrasts
fsl_surf_model.inputs.level1design.interscan_interval = exp.TR
fsl_surf_model.inputs.level1design.bases = exp.fsl_bases

selectsurfcontrast = fsl_surf_model.get_node("selectcontrast")
selectsurfcontrast.iterables = ('index', [[i] for i in range(len(contrasts))])


firstlevel.connect([
    (infosource, fsl_vol_model, 
        [(("subject_id", fuf.subjectinfo, exp), "modelspec.subject_info"),
         ("subject_id", "modelspec.subject_id")]),
    (infosource, fsl_surf_model, 
        [(("subject_id", fuf.subjectinfo, exp), "modelspec.subject_info"),
         ("subject_id", "modelspec.subject_id")]),
    (preproc, fsl_vol_model, 
        [("volhighpass.out_file", "modelspec.functional_runs"),
         ("art.outlier_files", "modelspec.outlier_files"),
         ("art.intensity_files", "stimcorr.intensity_values"),
         ("realign.par_file", "stimcorr.realignment_parameters"),
         ("realign.par_file", "modelspec.realignment_parameters"),
         ("volhighpass.out_file","modelestimate.in_file"),
         (("meanfunc2.out_file", lambda x: x[0]),"overlaystats.background_image"),
         (("meanfunc2.out_file", lambda x: x[0]),"overlayresidual.background_image"),
         ("reg2struct.out_reg_file", "xfmcopes.reg_file"),
         ("reg2struct.out_reg_file", "xfmvarcopes.reg_file"),
        ]),
    (preproc, fsl_surf_model, 
        [("surfhighpass.out_file", "modelspec.functional_runs"),
         ("art.outlier_files", "modelspec.outlier_files"),
         ("art.intensity_files", "stimcorr.intensity_values"),
         ("realign.par_file", "stimcorr.realignment_parameters"),
         ("realign.par_file", "modelspec.realignment_parameters"),
         ("surfhighpass.out_file","modelestimate.in_file"),
         (("meanfunc2.out_file", lambda x: x[0]),"overlaystats.background_image"),
         (("meanfunc2.out_file", lambda x: x[0]),"overlayresidual.background_image"),
         ("reg2struct.out_reg_file", "xfmcopes.reg_file"),
         ("reg2struct.out_reg_file", "xfmvarcopes.reg_file"),
        ]),
    ])


if args.runspm:
    spm_modelfit.inputs.modelspec.time_repetition = exp.TR
    spm_modelfit.inputs.modelspec.high_pass_filter_cutoff = exp.hpcutoff
    spm_modelfit.inputs.modelspec.input_units = exp.units
    spm_modelfit.inputs.modelspec.output_units = exp.units

    spm_modelfit.inputs.contrastestimate.contrasts = contrasts
    spm_modelfit.inputs.level1design.interscan_interval = exp.TR
    spm_modelfit.inputs.level1design.bases = exp.spm_bases
    spm_modelfit.inputs.level1design.timing_units = exp.units
    
    firstlevel.connect([(infosource, spm_modelfit, 
                            [(("subject_id", fuf.subjectinfo, exp), "modelspec.subject_info"),
                             ("subject_id", "modelspec.subject_id")]),
                        (preproc, spm_modelfit, 
                            [("unzip_intnorm.out_file", "modelspec.functional_runs"),
                             ("art.outlier_files", "modelspec.outlier_files"),
                             ("realign.par_file", "modelspec.realignment_parameters")])
                        ])

if exp.nruns > 1:
    firstlevel.connect(
        [(preproc, vol_fixed_fx, 
           [("dilatemask.out_file", "flameo.mask_file"),
           (("meanfunc2.out_file", lambda x: x[0]), "overlayflame.background_image")]),
         (fsl_vol_model, vol_fixed_fx,
           [(("xfmcopes.transformed_file", fuf.sort_copes),"copemerge.in_files"),
            (("xfmvarcopes.transformed_file", fuf.sort_copes),"varcopemerge.in_files"),
            (("contrastestimate.copes", lambda x: len(x)),"l2model.num_copes")]),
         (preproc, surf_fixed_fx, 
           [("dilatemask.out_file", "flameo.mask_file"),
           (("meanfunc2.out_file", lambda x: x[0]), "overlayflame.background_image")]),
         (fsl_surf_model, surf_fixed_fx,
           [(("xfmcopes.transformed_file", fuf.sort_copes),"copemerge.in_files"),
            (("xfmvarcopes.transformed_file", fuf.sort_copes),"varcopemerge.in_files"),
            (("contrastestimate.copes", lambda x: len(x)),"l2model.num_copes")]),
        ])


# File crashdumps by date
datestamp = str(datetime.now())[:10]
codepath = os.path.split(os.path.abspath(__file__))[0]
crashdir = os.path.abspath("%s/crashdumps/%s" % (codepath, datestamp))
if not os.path.isdir(crashdir):    
    os.makedirs(crashdir)
firstlevel.config = dict(crashdump_dir=crashdir) 

""" Setup the output """
working_output = os.path.join(
    os.path.abspath("../Analysis/NiPype/workingdir"), args.paradigm)
firstlevel.base_dir = working_output

output_base = os.path.join(os.path.abspath("../Analysis/NiPype"), args.paradigm)

datasink = pe.Node(interface=nio.DataSink(), 
                   name="datasink")
datasink.inputs.base_directory = output_base
datasink.overwrite = True

reportsub = []
imagesub = []

for r in range(exp.nruns):
    runstr = "run_%d"%(r+1)
    for type in ["rotations","translations"]:
        reportsub.append(("_plot_type_%s/_plotmotion%d"%(type,r),runstr))
    imagesub.append(("_realign%d"%r,""))
    reportsub.append(("_plotdisplacement%d"%r,runstr))
    reportsub.append(("_modelgen%d"%r,runstr))
    reportsub.append(("run%d"%r,"design"))
    imagesub.append(("_modelestimate%d"%r,runstr))
    reportsub.append(("_sliceresidual%d"%r,runstr))
    imagesub.append(("_highpass%d"%r,""))
for con, conparams in enumerate(contrasts):
    for r in range(exp.nruns):
        reportsub.append(("_index_%d/_slicestats%d/zstat%d_"%(con,r,con+1),
                          "run_%d/%s_"%(r+1,conparams[0])))
    reportsub.append(("_sliceflame%d/zstat1_"%con, "%s_"%conparams[0]))
    imagesub.append(("_flameo%d"%con,"%s"%conparams[0]))
reportsub.reverse()
imagesub.reverse()

mergesubs = pe.Node(interface=util.Merge(numinputs=5),
                       name="mergesubstitutes")

firstlevel.connect([(preproc, mergesubs,
                       [(("volhighpass.out_file",fuf.sub,"preproc-vol_func.nii.gz"),"in1"),
                        (("surfhighpass.out_file",fuf.sub,"preproc-surf_func.nii.gz"),"in2"),
                        (("art.outlier_files",fuf.sub,"outliers.txt"),"in3"),
                        (("realign.par_file",fuf.sub,"motion_params.par"),"in4"),
                        (("dilatemask.out_file",fuf.sub,"func_mask.nii.gz",False),"in5")])
                    ])

report = pe.Node(interface=nio.DataSink(),
                 name="report")
report.inputs.base_directory = os.path.join(output_base, "report")
report.inputs.substitutions = reportsub
report.inputs.parameterization = True
report.overwrite = True

firstlevel.connect([(infosource, datasink, 
                        [("subject_id", "container"),
                         (("subject_id", lambda x: "_subject_id_%s"%x), "strip_dir")]),
                    (mergesubs, datasink,
                        [(("out", lambda x: imagesub + x), "substitutions")]),
                    (preproc, datasink,
                        [("volhighpass.out_file", "preproc.@functional_runs"),
                         ("surfhighpass.out_file", "preproc.@functional_runs"),
                         ("art.outlier_files", "preproc.@outlier_files"),
                         ("realign.par_file", "preproc.@realign_parameters"),
                         ("dilatemask.out_file", "preproc.@func_mask")]),
                    (fsl_vol_model, datasink,
                        [("modelestimate.results_dir", "level1.model.@results"),
                         ("contrastestimate.tstats", "level1.contrasts.@T"),
                         ("contrastestimate.zstats", "level1.contrasts.@Z"),
                         ("contrastestimate.copes", "level1.contrasts.@copes"),
                         ("contrastestimate.varcopes", "level1.contrasts.@varcopes")]),
                    (fsl_surf_model, datasink,
                        [("modelestimate.results_dir", "level1.model.@results"),
                         ("contrastestimate.tstats", "level1.contrasts.@T"),
                         ("contrastestimate.zstats", "level1.contrasts.@Z"),
                         ("contrastestimate.copes", "level1.contrasts.@copes"),
                         ("contrastestimate.varcopes", "level1.contrasts.@varcopes")]),
                    (infosource, report,
                        [("subject_id", "container"),
                         (("subject_id", lambda x: "_subject_id_%s"%x), "strip_dir")]),
                    (preproc, report,
                        [("plotmotion.out_file", "preproc.@motionplots"),
                         ("plotdisplacement.out_file", "preproc.@displacementplots")]),
                    (fsl_vol_model, report,
                        [("modelgen.design_image", "level1.@design"),
                         ("modelgen.design_cov", "level1.@covariance"),
                         ("slicestats.out_file", "level1.@zstats"),
                         ("sliceresidual.out_file", "level1.@sigmasquared")]),
                    (fsl_surf_model, report,
                        [("modelgen.design_image", "level1.@design"),
                         ("modelgen.design_cov", "level1.@covariance"),
                         ("slicestats.out_file", "level1.@zstats"),
                         ("sliceresidual.out_file", "level1.@sigmasquared")]),
                    ])
if exp.nruns > 1:                    
    firstlevel.connect([(vol_fixed_fx, datasink,
                            [("flameo.stats_dir", "fixedfx.@zstats")]),
                        (vol_fixed_fx, report,
                            [("sliceflame.out_file", "fixedfx.@zstats")]),
                        (surf_fixed_fx, datasink,
                            [("flameo.stats_dir", "fixedfx.@zstats")]),
                        (surf_fixed_fx, report,
                            [("sliceflame.out_file", "fixedfx.@zstats")]),
                        ])
if args.runspm:
    firstlevel.connect([(spm_modelfit, datasink,
                            [("contrastestimate.con_images","SPM.contrasts.@con"),
                             ("contrastestimate.spmT_images","SPM.contrasts.@T")])
                        ])

# Run the script
if args.run:
    firstlevel.run(inseries=args.inseries)
    logdir = "/mindhive/gablab/fluid/NiPype_Code/log_archive/%s"%datestamp
    if not os.path.isdir(logdir):
        os.mkdir(logdir)
    timestamp = str(datetime.now())[11:16].replace(":","-")
    fulllog = open("%s/%s_%s.log"%(logdir,args.paradigm,timestamp),"w")
    for lf in ["pypeline.log%s"%n for n in [".4",".3",".2",".1",""]]:
        if os.path.isfile(lf):
            fulllog.write(open(lf).read())
            os.remove(lf)
    fulllog.close()
if args.write_graph:
    firstlevel.write_graph(graph2use="flat")
