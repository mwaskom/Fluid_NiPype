"""
    Nipype experiment module for the NBack paradigm.
"""
import re

template_args = dict(timeseries=[["subject_id",["IQ_run?"]]])

source_template = "%s/bold/%s.nii.gz"

hpcutoff = 128
TR = 2.
units = "secs"

nruns = 1

fsl_bases = {"dgamma":{"derivs":False}}
spm_bases = {"hrf":[0,0]}

parfile_base_dir = "/mindhive/gablab/fluid/Data"
parfile_template = "%s/parfiles/IQ_r%d_d1_%s_%s.txt"
parfile_args = ["subject_id", "run_number", "event", "subject_id"]

events = ["easy", "hard"]

cont01 = ["easy", "T", events, [1,0]]
cont02 = ["hard", "T", events, [0,1]]
cont03 = ["easy-hard", "T", events, [1,-1]]
cont04 = ["hard-easy", "T", events, [-1,1]]

convars = [var for var in dir() if re.match("cont\d+",var)]
convars.sort()

contrasts = [locals()[con] for con in convars]
