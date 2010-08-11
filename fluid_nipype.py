# vi: set ft=python sts=4 ts=4 sw=4 et:
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
"""
    Nipype interface script for the Fluid Intelligence project.
"""

import os
import sys
import argparse
from copy import deepcopy
from datetime import datetime

import nipype.pipeline.engine as pe
#import new_io as nio
import nipype.interfaces.io as nio
import nipype.interfaces.spm as spm
import nipype.interfaces.freesurfer as fs

# Catch when we"re not using the right environment
if not pe.__file__.startswith("/software/python/nipype0.3"):
    sys.exit("ERROR: Not using nipype0.3")

from fluid_preproc import preproc
from fluid_fsl_model import fsl_modelfit
from fluid_spm_model import spm_modelfit
from fluid_fsl_fixed_fx import fixed_fx

import fluid_utility_funcs as fuf

""" Handle command line arguments to control the analysis """
parser = argparse.ArgumentParser(description="Main interface for GFluid NiPype code.")

""" Define the paradigm.  We"ll eventually get this from the command line."""
paradigm = "mot_block"

""" Dynamically import the experiment file """
exp = __import__("%s_experiment" % paradigm)

""" Subjects.  This won"t stay hardcorded like this """
subject_list = ["SMARTER_SP17"] # , "SMARTER_SP15"]

""" Define the level 1 pipeline"""
firstlevel = pe.Workflow(name= "level1")


""" Tell infosource to iterate over all subjects """
exp.infosource.iterables = ("subject_id", subject_list)

""" Add experiment-sepcific information to the workflows """

preproc.inputs.highpass.suffix = "_hpf"
preproc.inputs.highpass.op_string = "-bptf %s -1" % (exp.hpcutoff/exp.TR)

preproc.inputs.smooth.fwhm = 6.

fsl_modelfit.inputs.modelspec.time_repetition = exp.TR
fsl_modelfit.inputs.modelspec.high_pass_filter_cutoff = exp.hpcutoff
fsl_modelfit.inputs.modelspec.input_units = exp.units
fsl_modelfit.inputs.modelspec.output_units = exp.units

spm_modelfit.inputs.modelspec.time_repetition = exp.TR
spm_modelfit.inputs.modelspec.high_pass_filter_cutoff = exp.hpcutoff
spm_modelfit.inputs.modelspec.input_units = exp.units
spm_modelfit.inputs.modelspec.output_units = exp.units

contrast_names = [n for n in exp.__dict__.keys() if "cont" in n] # Python magic
contrast_names.sort()
contrasts = []
for name in contrast_names:
    for k, v in exp.__dict__.items():
        if k == name:
            contrasts.append(v)

fsl_modelfit.inputs.level1design.contrasts = contrasts
fsl_modelfit.inputs.level1design.interscan_interval = exp.TR
fsl_modelfit.inputs.level1design.bases = exp.fsl_bases

spm_modelfit.inputs.contrastestimate.contrasts = contrasts
spm_modelfit.inputs.level1design.interscan_interval = exp.TR
spm_modelfit.inputs.level1design.bases = exp.spm_bases
spm_modelfit.inputs.level1design.timing_units = exp.units

""" Connect the workflows """
firstlevel.connect([(exp.infosource, exp.datasource, 
                        [("subject_id", "subject_id")]),
                    (exp.datasource, preproc, 
                        [("struct", "inputspec.struct"),
                         ("func", "inputspec.func")]),
                    (exp.infosource, fsl_modelfit, 
                        [(("subject_id", fuf.subjectinfo, exp), "modelspec.subject_info"),
                         ("subject_id", "modelspec.subject_id")]),
                    (exp.infosource, spm_modelfit, 
                        [(("subject_id", fuf.subjectinfo, exp), "modelspec.subject_info"),
                         ("subject_id", "modelspec.subject_id")]),
                    (preproc, fsl_modelfit, 
                        [("highpass.out_file", "modelspec.functional_runs"),
                         ("art.outlier_files", "modelspec.outlier_files"),
                         ("realign.par_file", "modelspec.realignment_parameters"),
                         ("highpass.out_file","modelestimate.in_file")]),
                    (preproc, spm_modelfit, 
                        [("unzip_intnorm.out_file", "modelspec.functional_runs"),
                         ("art.outlier_files", "modelspec.outlier_files"),
                         ("realign.par_file", "modelspec.realignment_parameters")]),
                    (preproc, fixed_fx, 
                        [("coregister.out_file", "flameo.mask_file")]),
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

substitutions = []
for i, name in enumerate(contrasts):
    substitutions.append(("_flameo%d"%i,contrasts[i][0]))
datasink.inputs.substitutions = substitutions    

firstlevel.connect([(exp.infosource, datasink, 
                        [("subject_id", "container"),
                         (("subject_id", lambda x: "_subject_id_%s"%x), "strip_dir")]),
                    (spm_modelfit, datasink,
                        [("contrastestimate.con_images","SPM.contrasts.@con"),
                         ("contrastestimate.spmT_images","SPM.contrasts.@T")]),
                    (fsl_modelfit, datasink,
                        [("contrastestimate.tstats", "FSL.level1.@T"),
                         ("contrastestimate.zstats", "FSL.level1.@Z"),
                         ("contrastestimate.copes", "FSL.level1.@copes"),
                         ("contrastestimate.varcopes", "FSL.level1.@varcopes")]),
                    (fixed_fx, datasink,
                        [("flameo.stats_dir", "FSL.fixedfx.@stats")])
                    ])

# Run the script
if __name__ == "__main__":
    if not "--norun" in sys.argv:
        firstlevel.run()
    firstlevel.write_graph(graph2use = "flat")
