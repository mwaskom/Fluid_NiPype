# vi: set ft=python sts=4 ts=4 sw=4 et:
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
"""
    Datasource module for fluid intelligence project pipelines.
"""

import os
import sys

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util

data_dir = "/mindhive/gablab/fluid/fmri_MOT_IQ_nback_pilot/data"
subject_list = ["SMARTER_SP14"]

infosource = pe.Node(interface=util.IdentityInterface(fields=["subj_id"]),
                     name="infosource")


datasource = pe.Node(interface=nio.DataGrabber(infields=["subj_id"],
                                               outfields=["func", "struct"]),
                     name="datasource")

datasource.inputs.base_directory = data_dir
datasource.inputs.template = "%s/%s/%s.nii.gz"
