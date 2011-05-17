"""
    Nipype experiment module for the jittered MOT paradigm.
"""
import re

template_args = dict(timeseries=[["subject_id", "bold", ["MOT_Jitter_run?"]]])

source_template = "%s/%s/%s.nii.gz"

hpcutoff = 128
TR = 2.
units = "secs"

nruns = 1

fsl_bases = {"dgamma":{"derivs":False}}

parfile_base_dir = "/mindhive/gablab/fluid/Data"
parfile_template = "%s/parfiles/MOT_Jitter_r%d_d%d_%s_%s.txt"
parfile_args = ["subject_id", "run_number", "day", "event", "subject_id"]

events = ["speed1", "speed2", "speed3", "speed4", "resp"]

cont01 = ["slow",         "T", events, [1,0,0,0,0]]
cont02 = ["normal",       "T", events, [0,1,0,0,0]]
cont03 = ["fast",         "T", events, [0,0,1,0,0]]
cont04 = ["vfast",        "T", events, [0,0,0,1,0]]
cont05 = ["normal-slow",  "T", events, [-1,1,0,0,0]]
cont06 = ["fast-slow",    "T", events, [-1,0,1,0,0]]
cont07 = ["fast-normal",  "T", events, [0,-1,1,0,0]]
cont08 = ["vfast-slow",   "T", events, [-1,0,0,1,0]]
cont09 = ["vfast-normal", "T", events, [0,-1,0,1,0]]
cont10 = ["vfast-fast",   "T", events, [0,0,-1,1,0]]
cont11 = ["all_speeds",   "T", events, [1./4,1./4,1./4,1./4,0]]

convars = [var for var in dir() if re.match("cont\d+",var)]
convars.sort()

contrasts = [locals()[con] for con in convars]
