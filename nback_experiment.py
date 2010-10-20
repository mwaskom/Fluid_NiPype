"""
    Nipype experiment module for the NBack paradigm.
"""
import re

template_args = dict(timeseries=[["subject_id", "bold", ["NBack_run?"]]])

source_template = "%s/%s/%s.nii.gz"

hpcutoff = 128
TR = 2.
units = "secs"

nruns = 4

fsl_bases = {"dgamma":{"derivs":False}}
spm_bases = {"hrf":[0,0]}

parfile_template = "%s/parfiles/NBack_r%d_d1_%s_%s.txt"

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

convars = [var for var in dir() if re.match("cont\d+",var)]
convars.sort()

contrasts = [locals()[con] for con in convars]

