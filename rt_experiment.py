# vi: set ft=python sts=4 ts=4 sw=4 et:
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
"""
    Nipype experiment module for the simply/choice rt paradigm.
"""

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util

data_dir = "/mindhive/gablab/fluid/fmri_MOT_IQ_nback_pilot/data"

infosource = pe.Node(interface=util.IdentityInterface(fields=["subject_id"]),
                     name="infosource")


datasource = pe.Node(interface=nio.DataGrabber(infields=["subject_id"],
                                               outfields=["func", "struct"]),
                     name="datasource")

datainfo = dict(func=[["subject_id", ["TestRT_run1", "TestRT_run2"]],
                struct=[["subject_id", "mprage"]])

datasource.inputs.base_directory = data_dir
datasource.inputs.template = "%s/nii/%s.nii.gz"

datasource.inputs.template_args = datainfo

hpcutoff = 128
TR = 2.
units = "secs"

nruns = 2

fsl_bases = {"dgamma":{"derivs":False}}
spm_bases = {"hrf":[0,0]}

parfile_template = "%s/parfiles/RT%s_%s_r%d.txt"

names = ["Simple", "Choice", "Inst"]

cont01 = ["simple", "T", names, [1,0,0]]
cont02 = ["choice", "T", names, [0,1,0]]
cont03 = ["simple-choice", "T", names, [1,-1,0]]
cont04 = ["choice-simple", "T", names, [-1,1,0]]
