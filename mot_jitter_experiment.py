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
                                               outfields=["func", "target", "struct"]),
                     name="datasource")

datainfo = dict(func=[["subject_id", "MOT_Jitter_run1"]],
                target=[["subject_id", "ep2d_t1w"]],
                struct=[["subject_id", "mprage"]])

datasource.inputs.base_directory = data_dir
datasource.inputs.template = "%s/nii/%s.nii.gz"

datasource.inputs.template_args = datainfo


exclude_subjects = ["SMARTER_SP%02d"%i for i in [2, 7, 8, 9]]

hpcutoff = 128
TR = 2.
units = "secs"

nruns = 1

fsl_bases = {"dgamma":{"derivs":False}}
spm_bases = {"hrf":[0,0]}

parfile_template = "%s/parfiles/MOT_Jitter_r%d_%s_%s.txt"

names = ["speed1", "speed2", "speed3", "speed4", "resp"]

cont01 = ["slow", "T", names, [1,0,0,0,0]]
cont02 = ["normal", "T", names, [0,1,0,0,0]]
cont03 = ["fast", "T", names, [0,0,1,0,0]]
cont04 = ["vfast", "T", names, [0,0,0,1,0]]
cont05 = ["slow-normal", "T", names, [1,-1,0,0,0]]
cont06 = ["slow-fast", "T", names, [1,0,-1,0,0]]
cont07 = ["slow-vfast", "T", names, [1,0,0,-1,0]]
cont08 = ["normal-slow", "T", names, [-1,1,0,0,0]]
cont09 = ["normal-fast", "T", names, [0,1,-1,0,0]]
cont10 = ["normal-vfast", "T", names, [0,1,0,-1,0]]
cont11 = ["fast-slow", "T", names, [-1,0,1,0,0]]
cont12 = ["fast-normal", "T", names, [0,-1,1,0,0]]
cont13 = ["fast-vfast", "T", names, [0,0,1,-1,0]]
cont14 = ["vfast-slow", "T", names, [-1,0,0,1,0]]
cont15 = ["vfast-normal", "T", names, [0,-1,0,1,0]]
cont16 = ["vfast-fast", "T", names, [0,0,-1,1,0]]

