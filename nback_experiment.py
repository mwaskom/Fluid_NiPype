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
                                               outfields=["func", "target", "struct"]),
                     name="datasource")

datainfo = dict(func=[["subject_id", "nii", 
                          ["NBack_run1"]]], #, "NBack_run2", "NBack_run3", "NBack_run4"]]],
                target=[["subject_id", "nii", "ep2d_t1w"]],
                struct=[["subject_id", "nii", "mprage"]])

datasource.inputs.base_directory = data_dir
datasource.inputs.template = "%s/%s/%s.nii.gz"

datasource.inputs.template_args = datainfo

hpcutoff = 128
TR = 2.
units = "secs"

nruns = 1

fsl_bases = {"dgamma":{"derivs":False}}
spm_bases = {"hrf":[0,0]}

parfile_template = "%s/parfiles/NBack_r%d_%s_%s.txt"

names = ["zero","easy","medium","hard","inst"]

cont01 = ["zero", "T", names, [1,0,0,0,0]]
cont02 = ["easy", "T", names, [0,1,0,0,0]]
cont03 = ["medium", "T", names, [0,0,1,0,0]]
cont04 = ["hard", "T", names, [0,0,0,1,0]]
cont05 = ["easy-zero", "T", names, [-1,1,0,0,0]]
cont06 = ["medium-zero", "T", names, [-1,0,1,0,0]]
cont07 = ["medium-easy", "T", names, [0,-1,1,0,0]]
cont08 = ["hard-zero", "T", names, [-1,0,0,1,0]]
cont09 = ["hard-easy", "T", names, [0,-1,0,1,0]]
cont10 = ["hard-medium", "T", names, [0,0,-1,1,0]]
cont11 = ["load-zero", "T", names, [-1,1./3,1./3,1./3,0]]
