"""
    Nipype experiment module for the jittered MOTT paradigm.
"""
import re

template_args = dict(func=[["subject_id", "bold", "MOT_Block_run1"]],
                     struct=[["subject_id", "structural", "mprage"]])

source_template = "%s/%s/%s.nii.gz"

hpcutoff = 128
TR = 2.
units = "secs"

nruns = 1

fsl_bases = {"dgamma":{"derivs":False}}
spm_bases = {"hrf":[0,0]}

parfile_template = "%s/parfiles/MOT_Block_r%d_d1_%s_%s.txt"

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

convars = [var for var in dir() if re.match("cont\d+",var)]
convars.sort()

contrasts = [locals()[con] for con in convars]
