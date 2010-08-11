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

datainfo = dict(func=[["subject_id", "nii", ["NBack_run1", "NBack_run2"]]],
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

parfile_template = "%s/parfiles/NBack_%s_d1_r%d_%s.txt"

names = ["0back","1back","2back","4back","inst"]

cont01 = ["0Back", "T", names, [1,0,0,0,0]]
cont02 = ["1Back", "T", names, [0,1,0,0,0]]
cont03 = ["2Back", "T", names, [0,0,1,0,0]]
cont04 = ["4Back", "T", names, [0,0,0,1,0]]
cont05 = ["1B-0B", "T", names, [-1,1,0,0,0]]
cont06 = ["2B-0B", "T", names, [-1,0,1,0,0]]
cont07 = ["2B-1B", "T", names, [0,-1,1,0,0]]
cont08 = ["4B-0B", "T", names, [-1,0,0,1,0]]
cont09 = ["4B-1B", "T", names, [0,-1,0,1,0]]
cont10 = ["4B-2B", "T", names, [0,0,-1,1,0]]
cont11 = ["421B-0B", "T", names, [-1,1./3,1./3,1./3,0]]
