"""
    Nipype experiment module for the NBack paradigm.
"""

template_args = dict(func=[["subject_id", "nii", 
                             ["NBack_run1", "NBack_run2", "NBack_run3", "NBack_run4"]]],
                     target=[["subject_id", "nii", "ep2d_t1w"]],
                     struct=[["subject_id", "nii", "mprage"]])

source_template = "%s/%s/%s.nii.gz"

hpcutoff = 128
TR = 2.
units = "secs"

nruns = 4

exclude_subjects = ["SMARTER_SP%02d"%i for i in range(17)]

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
