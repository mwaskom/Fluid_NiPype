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
if not pe.__file__.startswith("/software/python/nipype0.3"):
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
    subject_list = ["SMARTER_SP%02d"%(i+1) for i in range(23)]
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

contrast_names = [n for n in exp.__dict__.keys() if re.match("cont\d*", n)] # Python magic
contrast_names.sort()
contrasts = []
for name in contrast_names:
    for k, v in exp.__dict__.items():
        if k == name:
            contrasts.append(v)

fsl_modelfit.inputs.modelspec.time_repetition = exp.TR
fsl_modelfit.inputs.modelspec.high_pass_filter_cutoff = exp.hpcutoff
fsl_modelfit.inputs.modelspec.input_units = exp.units
fsl_modelfit.inputs.modelspec.output_units = exp.units

fsl_modelfit.inputs.level1design.contrasts = contrasts
fsl_modelfit.inputs.level1design.interscan_interval = exp.TR
fsl_modelfit.inputs.level1design.bases = exp.fsl_bases

firstlevel.connect([(infosource, fsl_modelfit, 
                        [(("subject_id", fuf.subjectinfo, exp), "modelspec.subject_info"),
                         ("subject_id", "modelspec.subject_id")]),
                    (preproc, fsl_modelfit, 
                        [("highpass.out_file", "modelspec.functional_runs"),
                         ("art.outlier_files", "modelspec.outlier_files"),
                         ("realign.par_file", "modelspec.realignment_parameters"),
                         ("highpass.out_file","modelestimate.in_file")]),
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
                       [("binarize_func.out_file", "flameo.mask_file")]),
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
    os.path.abspath("../nipype_output/pilots/workingdir"), paradigm)
firstlevel.base_dir = working_output

datasink = pe.Node(interface=nio.DataSink(), 
                   name="datasink")
datasink.inputs.base_directory = os.path.join(
    os.path.abspath("../nipype_output/pilots"), paradigm)

if exp.nruns == 1:
    datasink.inputs.parameterization = False

substitutions = []
for i, name in enumerate(contrasts):
    substitutions.append(("_flameo%d"%i,"%s"%contrasts[i][0]))
for i in range(exp.nruns):
    substitutions.append(("_modelestimate%s"%i,"run_%s"%(i+1)))
substitutions.reverse()
datasink.inputs.substitutions = substitutions    

firstlevel.connect([(infosource, datasink, 
                        [("subject_id", "container"),
                         (("subject_id", lambda x: "_subject_id_%s"%x), "strip_dir")]),
                    (fsl_modelfit, datasink,
                        [("contrastestimate.tstats", "FSL.level1.@T"),
                         ("contrastestimate.zstats", "FSL.level1.@Z"),
                         ("contrastestimate.copes", "FSL.level1.@copes"),
                         ("contrastestimate.varcopes", "FSL.level1.@varcopes")])
                    ])
if exp.nruns > 1:                    
    firstlevel.connect([(fixed_fx, datasink,
                            [("flameo.stats_dir", "FSL.fixedfx.@stats")])])
if do_spm:
    firstlevel.connect([(spm_modelfit, datasink,
                            [("contrastestimate.con_images","SPM.contrasts.@con"),
                             ("contrastestimate.spmT_images","SPM.contrasts.@T")])
                        ])

# Run the script
if __name__ == "__main__":
    if not "--norun" in sys.argv:
        firstlevel.run()
    if not "--nograph" in sys.argv:
        firstlevel.write_graph(graph2use = "flat")
