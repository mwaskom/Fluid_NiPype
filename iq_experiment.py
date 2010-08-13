# vi: set ft=python sts=4 ts=4 sw=4 et:
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
"""
    Nipype experiment module for the NBack paradigm.
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

datainfo = dict(func=[["subject_id", "nii", "IQ"]],
            struct=[["subject_id", "nii", "mprage"]])

datasource.inputs.base_directory = data_dir
datasource.inputs.template = "%s/%s/%s.nii.gz"

datasource.inputs.template_args = datainfo

hpcutoff = 128
TR = 2.
units = "secs"

nruns = 2

fsl_bases = {"dgamma":{"derivs":False}}
spm_bases = {"hrf":[0,0]}

parfile_template = "%s/parfiles/IQ_r%d_%s_%s.txt"

names = ["Easy", "Hard"]

cont01 = ["easy", "T", names, [1,0]]
cont02 = ["hard", "T", names, [0,1]]
cont03 = ["easy-hard", "T", names, [1,-1]]
cont04 = ["hard-easy", "T", names, [-1,1]]
