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
preproc.inputs.fssource.subjects_dir = "/mindhive/gablab/fluid/Data"

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
                        [("subject_id", "inputspec.subject_id"),
                         ("subject_id", "artspec.subject_id"),
                         (("subject_id", fuf.subjectinfo, exp), "artspec.subject_info")]),
                    (datasource, preproc, 
                        [("func", "inputspec.func")])
                    ])

contrasts = exp.contrasts

preproc.inputs.artspec.time_repetition = exp.TR
preproc.inputs.artspec.high_pass_filter_cutoff = exp.hpcutoff
preproc.inputs.artspec.input_units = exp.units
preproc.inputs.artspec.output_units = exp.units

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
    (preproc, fsl_vol_model, 
        [("volhighpass.out_file", "modelspec.functional_runs"),
         ("art.outlier_files", "modelspec.outlier_files"),
         ("realign.par_file", "modelspec.realignment_parameters"),
         ("volhighpass.out_file","modelestimate.in_file"),
         (("meanfunc2.out_file", lambda x: x[0]),"overlaystats.background_image"),
         (("meanfunc2.out_file", lambda x: x[0]),"overlayresidual.background_image"),
         ("reg2struct.out_reg_file", "xfmcopes.reg_file"),
         ("reg2struct.out_reg_file", "xfmvarcopes.reg_file"),
        ]),
    (infosource, fsl_surf_model, 
        [(("subject_id", fuf.subjectinfo, exp), "modelspec.subject_info"),
         ("subject_id", "modelspec.subject_id")]),
    (preproc, fsl_surf_model, 
        [("surfhighpass.out_file", "modelspec.functional_runs"),
         ("art.outlier_files", "modelspec.outlier_files"),
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
           [("mask2nii.out_file", "flameo.mask_file"),
            ("meanmean.concatenated_file", "overlayflame.background_image")]),
         (fsl_vol_model, vol_fixed_fx,
           [("xfmcopes.transformed_file","copemerge.in_files"),
            ("xfmvarcopes.transformed_file", "varcopemerge.in_files"),
            (("xfmcopes.transformed_file", lambda x: len(x)), "l2model.num_copes")]),
         (preproc, surf_fixed_fx, 
           [("mask2nii.out_file", "flameo.mask_file"),
           ("meanmean.concatenated_file", "overlayflame.background_image")]),
         (fsl_surf_model, surf_fixed_fx,
           [("xfmcopes.transformed_file","copemerge.in_files"),
            ("xfmvarcopes.transformed_file","varcopemerge.in_files"),
            (("xfmcopes.transformed_file", lambda x: len(x)),"l2model.num_copes")]),
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

l1reportsub = []
l2reportsub = []
sinksub = []

for r in range(exp.nruns):
    runstr = "run_%d"%(r+1)
    for type in ["rotations","translations"]:
        l1reportsub.append(("_plot_type_%s/_plotmotion%d"%(type,r),runstr))
    sinksub.append(("_realign%d"%r,""))
    l1reportsub.append(("_plotdisplacement%d"%r,runstr))
    l1reportsub.append(("_modelgen%d"%r,runstr))
    l1reportsub.append(("run%d"%r,"design"))
    sinksub.append(("_modelestimate%d"%r,runstr))
    l1reportsub.append(("_sliceresidual%d"%r,runstr))
    sinksub.append(("_highpass%d"%r,""))
    sinksub.append(("_art%d"%r,""))
    sinksub.append(("_dilatemask%d"%r,""))
    sinksub.append(("_volhighpass%d"%r,""))
    sinksub.append(("_surfhighpass%d"%r,""))
for con, conparams in enumerate(contrasts):
    for r in range(exp.nruns):
        l1reportsub.append(("_index_%d/_slicestats%d/zstat%d_"%(con,r,con+1),
                            "run_%d/%s_"%(r+1,conparams[0])))
    l2reportsub.append(("_sliceflame%d/zstat1_"%con, "%s_"%conparams[0]))
    sinksub.append(("_flameo%d"%con,"%s"%conparams[0]))
l1reportsub.reverse()
l2reportsub.reverse()
sinksub.reverse()

mergesubs = pe.Node(interface=util.Merge(numinputs=5),
                       name="mergesubstitutes")

firstlevel.connect([(preproc, mergesubs,
                       [(("volhighpass.out_file",fuf.sub,"preproc-vol_func.nii.gz"),"in1"),
                        (("surfhighpass.out_file",fuf.sub,"preproc-surf_func.nii.gz"),"in2"),
                        (("art.outlier_files",fuf.sub,"outliers.txt"),"in3"),
                        (("realign.par_file",fuf.sub,"motion_params.par"),"in4"),
                        (("dilatemask.out_file",fuf.sub,"func_mask.nii.gz"),"in5")])
                    ])

reportdir = os.path.join(output_base, "report")

volreport = pe.Node(interface=nio.DataSink(base_directory=reportdir),
                    overwrite=True,
                    name="volreport")

surfreport = pe.Node(interface=nio.DataSink(base_directory=reportdir),
                     overwrite=True,
                     name="surfreport")

firstlevel.connect([(infosource, datasink, 
                        [("subject_id", "container"),
                         (("subject_id", lambda x: "_subject_id_%s"%x), "strip_dir")]),
                    #(mergesubs, datasink,
                        #[(("out", lambda x: sinksub + x), "substitutions")]),
                    (preproc, datasink,
                        [("volhighpass.out_file", "preproc.@vol_runs"),
                         ("surfhighpass.out_file", "preproc.@surf_runs"),
                         ("art.outlier_files", "preproc.@outlier_files"),
                         ("realign.par_file", "preproc.@realign_parameters"),
                         ("dilatemask.out_file", "preproc.@func_mask")]),
                    (fsl_vol_model, datasink,
                        [("modelestimate.results_dir", "level1.volume.model.@results"),
                         ("contrastestimate.tstats", "level1.volume.contrasts.@T"),
                         ("contrastestimate.zstats", "level1.volume.contrasts.@Z"),
                         ("contrastestimate.copes", "level1.volume.contrasts.@copes"),
                         ("contrastestimate.varcopes", "level1.volume.contrasts.@varcopes")]),
                    (fsl_surf_model, datasink,
                        [("modelestimate.results_dir", "level1.surface.model.@results"),
                         ("contrastestimate.tstats", "level1.surface.contrasts.@T"),
                         ("contrastestimate.zstats", "level1.surface.contrasts.@Z"),
                         ("contrastestimate.copes", "level1.surface.contrasts.@copes"),
                         ("contrastestimate.varcopes", "level1.surface.contrasts.@varcopes")]),
                    (infosource, volreport,
                        [("subject_id", "container"),
                         (("subject_id", lambda x: "_subject_id_%s"%x), "strip_dir")]),
                    (infosource, volreport,
                        [("subject_id", "container"),
                         (("subject_id", lambda x: "_subject_id_%s"%x), "strip_dir")]),
                    (preproc, volreport,
                        [("plotmotion.out_file", "preproc.@motionplots"),
                         ("plotdisplacement.out_file", "preproc.@displacementplots")]),
                    (fsl_vol_model, volreport,
                        [("modelgen.design_image", "level1.@design"),
                         ("modelgen.design_cov", "level1.@covariance"),
                         ("slicestats.out_file", "level1.volume.@zstats"),
                         ("sliceresidual.out_file", "level1.volume.@sigmasquared")]),
                    (fsl_surf_model, surfreport,
                        [("slicestats.out_file", "level1.surface.@zstats"),
                         ("sliceresidual.out_file", "level1.surface.@sigmasquared")]),
                    ])
if exp.nruns > 1:                    
    firstlevel.connect([(infosource, volreport, 
                            [("subject_id", "container"),
                             (("subject_id", lambda x: "_subject_id_%s"%x), "strip_dir")]),
                        (infosource, surfreport, 
                            [("subject_id", "container"),
                             (("subject_id", lambda x: "_subject_id_%s"%x), "strip_dir")]),
                        (vol_fixed_fx, datasink,
                            [("flameo.stats_dir", "fixedfx.volume.@stats")]),
                        (vol_fixed_fx, volreport,
                            [("sliceflame.out_file", "fixedfx.volume.@report")]),
                        (surf_fixed_fx, datasink,
                            [("flameo.stats_dir", "fixedfx.surface.@stats")]),
                        (surf_fixed_fx, surfreport,
                            [("sliceflame.out_file", "fixedfx.surface.@report")]),
                        ])
if args.runspm and 0:
    firstlevel.connect([(spm_modelfit, datasink,
                            [("contrastestimate.con_images","SPM.contrasts.@con"),
                             ("contrastestimate.spmT_images","SPM.contrasts.@T")])
                        ])

# Run the script
if args.run:
    firstlevel.run(inseries=args.inseries)
    # Archive the logfile
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
