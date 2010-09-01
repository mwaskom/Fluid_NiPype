"""
    Experiment module for reaction-time paradigm
"""
import re

template_args = dict(func=[["subject_id", "bold", ["RT_run1", "RT_run2"]]],
                     target=[["subject_id", "structural", "ep2d_t1w"]],
                    struct=[["subject_id", "structrual", "mprage"]])

source_template = "%s/nii/%s.nii.gz"

hpcutoff = 128
TR = 2.
units = "secs"

nruns = 2

fsl_bases = {"dgamma":{"derivs":False}}
spm_bases = {"hrf":[0,0]}

parfile_template = "%s/parfiles/RT_r%d_%s_%s.txt"

names = ["simple", "choice", "inst"]

cont01 = ["simple", "T", names, [1,0,0]]
cont02 = ["choice", "T", names, [0,1,0]]
cont03 = ["simple-choice", "T", names, [1,-1,0]]
cont04 = ["choice-simple", "T", names, [-1,1,0]]

convars = [var for var in dir() if re.match("cont\d+",var)]
convars.sort()

contrasts = [locals()[con] for con in convars]
