# vi: set ft=python sts=4 ts=4 sw=4 et:
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
"""
    Nipype interface script for the Fluid Intelligence project.
"""

import os
import re
import sys
from copy import deepcopy
from datetime import datetime

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.spm as spm
import nipype.interfaces.freesurfer as fs

# Catch when we"re not using the right environment
if (not pe.__file__.startswith("/software/python/nipype0.3") 
    and not pe.__file__.startswith("/u2/mwaskom/nipype")):
    sys.exit("ERROR: Not using nipype0.3")

import argparse

from fluid_preproc import preproc
from fluid_fsl_model import fsl_modelfit
from fluid_spm_model import spm_modelfit
from fluid_fsl_fixed_fx import fixed_fx
from fluid_source import infosource, datasource

import fluid_utility_funcs as fuf

""" Handle command line arguments to control the analysis """
parser = argparse.ArgumentParser(description="Main interface for GFluid NiPype code.")

if "--spm" in sys.argv:
    do_spm = True
else:
    do_spm = False


""" Define the paradigm.  We"ll eventually get this from the command line."""
if "--nback" in sys.argv:
    paradigm = "nback"
elif "--motjitter" in sys.argv:
    paradigm = "mot_jitter"
elif "--motblock" in sys.argv:
    paradigm = "mot_block"
elif "--iq" in sys.argv:
    paradigm = "iq"
elif "--rt2" in sys.argv:
    paradigm = "rt_tworun"
elif "--rt" in sys.argv:
    paradigm = "rt"

""" Dynamically import the experiment file """
exp = __import__("%s_experiment" % paradigm)

""" Subjects.  This won"t stay hardcorded like this """
if hasattr(exp, "subject_list"):
    subject_list = exp.subject_list
else:
    subject_list = ["gf%02d"%id for id in [05, 14]]
if hasattr(exp, "exclude_subjects"):    
    subject_list = [subj for subj in subject_list if subj not in exp.exclude_subjects]

if "--subject" in sys.argv:
    subject_list = [sys.argv[sys.argv.index("--subject") + 1]]

""" Set experimental source attributes. """
datasource.inputs.template = exp.source_template
datasource.inputs.template_args = exp.template_args

exp.data_dir = datasource.inputs.base_directory

""" Define the level 1 pipeline"""
firstlevel = pe.Workflow(name= "level1")

""" Tell infosource to iterate over all subjects """
infosource.iterables = ("subject_id", subject_list)

""" Add experiment-sepcific information to the workflows """
preproc.inputs.highpass.suffix = "_hpf"
preproc.inputs.highpass.op_string = "-bptf %s -1" % (exp.hpcutoff/exp.TR)

preproc.inputs.smooth.fwhm = 6.

""" Connect the workflows """
firstlevel.connect([(infosource, datasource, 
                        [("subject_id", "subject_id")]),
                    (datasource, preproc, 
                        [("target", "inputspec.target"),
                         ("func", "inputspec.func")])
                   ])

contrasts = exp.contrasts

fsl_modelfit.inputs.modelspec.time_repetition = exp.TR
fsl_modelfit.inputs.modelspec.high_pass_filter_cutoff = exp.hpcutoff
fsl_modelfit.inputs.modelspec.input_units = exp.units
fsl_modelfit.inputs.modelspec.output_units = exp.units

fsl_modelfit.inputs.level1design.contrasts = contrasts
fsl_modelfit.inputs.level1design.interscan_interval = exp.TR
fsl_modelfit.inputs.level1design.bases = exp.fsl_bases

selectcontrast = fsl_modelfit.get_node("selectcontrast")
selectcontrast.iterables = ('index', [[i] for i in range(len(contrasts))])

firstlevel.connect([(infosource, fsl_modelfit, 
                        [(("subject_id", fuf.subjectinfo, exp), "modelspec.subject_info"),
                         ("subject_id", "modelspec.subject_id")]),
                    (preproc, fsl_modelfit, 
                        [("highpass.out_file", "modelspec.functional_runs"),
                         ("art.outlier_files", "modelspec.outlier_files"),
                         ("realign.par_file", "modelspec.realignment_parameters"),
                         ("highpass.out_file","modelestimate.in_file"),
                         (("meanfunc3.out_file", lambda x: x[0]),"overlaystats.background_image"),
                         (("meanfunc3.out_file", lambda x: x[0]),"overlayresidual.background_image")]),
                    ])


if do_spm:
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
                    [(preproc, fixed_fx, 
                       [("binarize_func.out_file", "flameo.mask_file"),
                       (("meanfunc3.out_file", lambda x: x[0]), "overlayflame.background_image")]),
                     (fsl_modelfit, fixed_fx,
                       [(("contrastestimate.copes", fuf.sort_copes),"copemerge.in_files"),
                        (("contrastestimate.varcopes", fuf.sort_copes),"varcopemerge.in_files"),
                        (("contrastestimate.copes", lambda x: len(x)),"l2model.num_copes")])
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
    os.path.abspath("../Analysis/NiPype/workingdir"), paradigm)
firstlevel.base_dir = working_output

output_base = os.path.join(os.path.abspath("../Analysis/NiPype"), paradigm)

datasink = pe.Node(interface=nio.DataSink(), 
                   name="datasink")
datasink.inputs.base_directory = output_base

reportsub = []
imagesub = []

for r in range(exp.nruns):
    reportsub.append(("_sliceresidual%d"%r,"run_%d"%(r+1)))
    imagesub.append(("_modelestimate%d"%r,"run_%d"%(r+1)))
for con, conparams in enumerate(contrasts):
    for r in range(exp.nruns):
        reportsub.append(("_index_%d/_slicestats%d/zstat%d_"%(con,r,con+1),
                          "run_%d/%s_"%(r+1,conparams[0])))
    reportsub.append(("_sliceflame%d/zstat1_"%con, "%s_"%conparams[0]))
    imagesub.append(("_flameo%d"%con,"%s"%conparams[0]))
reportsub.reverse()
imagesub.reverse()

datasink.inputs.substitutions = imagesub

report = pe.Node(interface=nio.DataSink(),
                 name="report")
report.inputs.base_directory = os.path.join(output_base, "report")
report.inputs.substitutions = reportsub
report.inputs.parameterization = True

firstlevel.connect([(infosource, datasink, 
                        [("subject_id", "container"),
                         (("subject_id", lambda x: "_subject_id_%s"%x), "strip_dir")]),
                    (fsl_modelfit, datasink,
                        [("modelestimate.results_dir", "level1.model.@results"),
                         ("contrastestimate.tstats", "level1.contrasts.@T"),
                         ("contrastestimate.zstats", "level1.contrasts.@Z"),
                         ("contrastestimate.copes", "level1.contrasts.@copes"),
                         ("contrastestimate.varcopes", "level1.contrasts.@varcopes")]),
                    (infosource, report,
                        [("subject_id", "container"),
                         (("subject_id", lambda x: "_subject_id_%s"%x), "strip_dir")]),
                    (fsl_modelfit, report,
                        [("slicestats.out_file", "level1.@zstats"),
                         ("sliceresidual.out_file", "level1.@sigmasquared")]),
                    ])
if exp.nruns > 1:                    
    firstlevel.connect([(fixed_fx, datasink,
                            [("flameo.stats_dir", "fixedfx.@zstats")]),
                        (fixed_fx, report,
                            [("sliceflame.out_file", "fixedfx.@zstats")])
                        ])
if do_spm:
    firstlevel.connect([(spm_modelfit, datasink,
                            [("contrastestimate.con_images","SPM.contrasts.@con"),
                             ("contrastestimate.spmT_images","SPM.contrasts.@T")])
                        ])

# Run the script
if __name__ == "__main__":
    if "--inseries" in sys.argv:
        inseries = True
    else:
        inseries = False
    if not "--norun" in sys.argv:
        firstlevel.run(inseries=inseries)
    if not "--nograph" in sys.argv:
        firstlevel.write_graph(graph2use = "flat")
