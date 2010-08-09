# vi: set ft=python sts=4 ts=4 sw=4 et:
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
"""
    Nipype interface script for the Fluid Intelligence project.
"""


import os
import sys
from datetime import datetime

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.spm as spm
import nipype.interfaces.freesurfer as fs

# Catch when we"re not using the right environment
if not pe.__file__.startswith("/software/python/nipype0.3"):
    sys.exit("ERROR: Not using nipype0.3")

from fluid_preproc import preproc
from fluid_datasource import data_dir, infosource, datasource
from fluid_modelfit import modelfit
from fluid_fixedfx import fixed_fx
from spm_tutorial2 import l1analysis as spm_analysis

# Subjects.  This won"t stay hardcorded like this
subject_list = ["SMARTER_SP15"]

# We"ll eventually get this from the command line
paradigm = "nback"

# Dynamically import the experiment file 
exp = __import__("%s_experiment" % paradigm)

# Get a reference to the highpass filter node
highpass = preproc.get_node("highpass")

# Set the highpass filter for the paradigm
highpass.inputs.suffix = "_hpf"
highpass.inputs.op_string = "-bptf %s -1" % (exp.hpcutoff/exp.TR)

# Get a reference to the smooth node
smooth = preproc.get_node("smooth")

# Set smoothing for the paradigm
smooth.inputs.fwhm = 6.

# Give the datasource node the template args for this paradigm
datasource.inputs.template_args = exp.datainfo

# Tell infosource to iterate over all subjects
infosource.iterables = ("subject_id", subject_list)

# Connect preproc and model flows
def sort_copes(files):
    numelements = len(files[0])
    outfiles = []
    for i in range(numelements):
        outfiles.insert(i,[])
        for j, elements in enumerate(files):
            outfiles[i].append(elements[i])
    return outfiles

def num_copes(files):
    return len(files)

# Convert the highpass filter output to .nii so SPM can read it
convert = pe.MapNode(interface=fs.MRIConvert(out_type="nii"), 
                     iterfield = "in_file",
                     name="convert")
preproc.connect(highpass, "out_file", convert, "in_file")

firstlevel = pe.Workflow(name="firstlevel")
firstlevel.connect([(preproc, modelfit, [("convert.out_file", "modelspec.functional_runs"),
                                         ("art.outlier_files", "modelspec.outlier_files"),
                                         ("realign.par_file", "modelspec.realignment_parameters"),
                                         ("highpass.out_file","modelestimate.in_file")]),
                    (preproc, fixed_fx, [("coregister.out_file", "flameo.mask_file")]),
                    (modelfit, fixed_fx,[(("conestimate.copes", sort_copes),"copemerge.in_files"),
                                         (("conestimate.varcopes", sort_copes),"varcopemerge.in_files"),
                                         (("conestimate.copes", num_copes),"l2model.num_copes"),
                                         ])
                    ])


# Subjectinfo function returns subj specific info

def parse_par_file(base_dir, subject_id, run_number, name, template):
    parfile = os.path.join(base_dir, template % (subject_id, subject_id, run_number, name))
    fid = open(parfile,"r")
    onsets = []
    durations = []
    for line in fid:
        line = line.split()
        onsets.append(float(line[0]))
        durations.append(float(line[1]))
    return onsets, durations

from nipype.interfaces.base import Bunch
from copy import deepcopy
def subjectinfo(subject_id):
    print "Subject ID: %s\n"%str(subject_id)
    output = []
    names = exp.names
    for r in range(exp.nruns):
        onsets = []
        durations = []
        for name in names:
            o, d = parse_par_file(data_dir, subject_id, r+1, name, exp.partemp)
            onsets.append(o)
            durations.append(d)
        output.insert(r,
                      Bunch(conditions=names,
                            onsets=deepcopy(onsets),
                            durations=deepcopy(durations),
                            amplitudes=None,
                            tmod=None,
                            pmod=None,
                            regressor_names=None,
                            regressors=None))
    return output


firstlevel.inputs.modelfit.modelspec.time_repetition = exp.TR
firstlevel.inputs.modelfit.modelspec.high_pass_filter_cutoff = exp.hpcutoff
firstlevel.inputs.modelfit.modelspec.input_units = "secs"
firstlevel.inputs.modelfit.modelspec.output_units = "secs"

contrast_names = [n for n in exp.__dict__.keys() if "cont" in n] # Python magic
contrast_names.sort()
contrasts = []
for name in contrast_names:
    for k, v in exp.__dict__.items():
        if k == name:
            contrasts.append(v)

firstlevel.inputs.modelfit.level1design.contrasts = contrasts
firstlevel.inputs.modelfit.level1design.interscan_interval = exp.TR
firstlevel.inputs.modelfit.level1design.bases = exp.bases

l1pipeline = pe.Workflow(name= "level1")


l1pipeline.connect([(infosource, datasource, [("subject_id", "subject_id")]),
                    (infosource, firstlevel, 
                        [(("subject_id", subjectinfo), "modelfit.modelspec.subject_info"),
                         ("subject_id", "modelfit.modelspec.subject_id")]),
                    (datasource, firstlevel, [("struct", "preproc.inputspec.struct"),
                                              ("func", "preproc.inputspec.func"),
                                              ]),
                    ])

# Add in the SPM analysis
"""Use :class:`nipype.interfaces.spm.EstimateModel` to determine the
parameters of the model.
"""
spm_level1 = pe.Workflow(name= "spm_level1")

spm_level1design = pe.Node(interface=spm.Level1Design(), name= "spm_level1design")
spm_level1design.inputs.bases              = {'hrf':{'derivs': [0,0]}}

spm_level1estimate = pe.Node(interface=spm.EstimateModel(), name="spm_level1estimate")
spm_level1estimate.inputs.estimation_method = {'Classical' : 1}

"""Use :class:`nipype.interfaces.spm.EstimateContrast` to estimate the
first level contrasts specified in a few steps above.
"""

spm_contrastestimate = pe.Node(interface = spm.EstimateContrast(), name="spm_contrastestimate")


spm_level1.connect([(spm_level1design, spm_level1estimate, [("spm_mat_file", "spm_mat_file")]),
                    (spm_level1estimate, spm_contrastestimate,[('spm_mat_file','spm_mat_file'),
                                                               ('beta_images','beta_images'),
                                                               ('residual_image','residual_image')]),
                  ])
spm_level1.inputs.spm_level1design.timing_units = "secs"
spm_level1.inputs.spm_level1design.interscan_interval = exp.TR
spm_level1.inputs.spm_contrastestimate.contrasts = contrasts


l1pipeline.connect(modelfit, "modelspec.session_info", spm_level1, "spm_level1design.session_info")


# File crashdumps by date
datestamp = str(datetime.now())[:10]
codepath = os.path.split(os.path.abspath(__file__))[0]
crashdir = os.path.abspath("%s/crashdumps/%s" % (codepath, datestamp))
if not os.path.isdir(crashdir):    
    os.makedirs(crashdir)
l1pipeline.config = dict(crashdump_dir=crashdir) 

# Setup the preproc working directory
working_output = os.path.abspath("../nipype_output/workingdir/singlerun/run2")
l1pipeline.base_dir = working_output

# Datasink
# --------

datasink = pe.Node(interface=nio.DataSink(), name="datasink")
datasink.inputs.base_directory = os.path.abspath("../nipype_output/")

l1pipeline.connect([(infosource, datasink,[("subject_id", "container")]),
                    (spm_level1, datasink,
                        [("spm_contrastestimate.con_images","SPM.contrasts.@con"),
                         ("spm_contrastestimate.spmT_images","SPM.contrasts.@T")]),
                    (firstlevel, datasink,
                        [("modelfit.conestimate.tstats", "FSL.level1.@T"),
                         ("modelfit.conestimate.zstats", "FSL.level1.@Z"),
                         ("modelfit.conestimate.copes", "FSL.level1.@copes"),
                         ("modelfit.conestimate.varcopes", "FSL.level1.@varcopes"),
                         ("fixedfx.flameo.zstats", "FSL.fixedfx.@Z"),
                         ("fixedfx.flameo.tstats", "FSL.fixedfx.@T"),
                         ("fixedfx.flameo.copes", "FSL.fixedfx.@copes"),
                         ("fixedfx.flameo.var_copes", "FSL.fixedfx.@varcopes")])])


# Run the script
if __name__ == "__main__":
    if not "--norun" in sys.argv:
        l1pipeline.run()
    l1pipeline.write_graph(graph2use = "flat")
