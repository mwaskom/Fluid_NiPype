"""
    Nipype experiment module for the NBack paradigm.
"""
import re

template_args = dict(func=[["subject_id", "IQ_run1"]],
                    target=[["subject_id", "ep2d_t1w"]],
                    struct=[["subject_id", "mprage"]])

source_template = "%s/nii/%s.nii.gz"

hpcutoff = 128
TR = 2.
units = "secs"

nruns = 1

fsl_bases = {"dgamma":{"derivs":False}}
spm_bases = {"hrf":[0,0]}

parfile_template = "%s/parfiles/IQ_r%d_d1_%s_%s.txt"

names = ["easy", "hard"]

cont01 = ["easy", "T", names, [1,0]]
cont02 = ["hard", "T", names, [0,1]]
cont03 = ["easy-hard", "T", names, [1,-1]]
cont04 = ["hard-easy", "T", names, [-1,1]]

convars = [var for var in dir() if re.match("cont\d+",var)]
convars.sort()

contrasts = [globals()[con] for con in convars]