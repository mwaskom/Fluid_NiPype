# vi: set ft=python sts=4 ts=4 sw=4 et:
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
"""
    Nipype experiment module for the jittered MOTT paradigm.
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

datainfo = dict(func=[["subject_id", "MOT_Block_run1"]],
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

parfile_template = "%s/parfiles/MOT%s_Bitter_r%d_%s.txt"

names = ["speed1", "speed2", "speed3", "speed4", "resp"]

cont01 = ["S1", "T", names, [1,0,0,0,0]]
cont02 = ["S2", "T", names, [0,1,0,0,0]]
cont03 = ["S3", "T", names, [0,0,1,0,0]]
cont04 = ["S4", "T", names, [0,0,0,1,0]]
cont05 = ["S1-S2", "T", names, [1,-1,0,0,0]]
cont06 = ["S1-S3", "T", names, [1,0,-1,0,0]]
cont07 = ["S1-S4", "T", names, [1,0,0,-1,0]]
cont08 = ["S2-S1", "T", names, [-1,1,0,0,0]]
cont09 = ["S2-S3", "T", names, [0,1,-1,0,0]]
cont10 = ["S2-S4", "T", names, [0,1,0,-1,0]]
cont11 = ["S3-S1", "T", names, [-1,0,1,0,0]]
cont12 = ["S3-S2", "T", names, [0,-1,1,0,0]]
cont13 = ["S3-S4", "T", names, [0,0,1,-1,0]]
cont14 = ["S4-S1", "T", names, [-1,0,0,1,0]]
cont15 = ["S4-S2", "T", names, [0,-1,0,1,0]]
cont16 = ["S4-S3", "T", names, [0,0,-1,1,0]]

