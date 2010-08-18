"""
    Nipype experiment module for the simple/choice reaction-time paradigm.
"""
template_args = dict(func=[["subject_id", "TestRT_run1"]],
                     target=[["subject_id", "ep2d_t1w"]],
                     struct=[["subject_id", "mprage"]])

source_template = "%s/nii/%s.nii.gz"

hpcutoff = 128
TR = 2.
units = "secs"

nruns = 1

exclude_subjects = ["SMARTER_SP%02d"%i for i in [1,2,3,7,8,9,10,11]]

fsl_bases = {"dgamma":{"derivs":False}}
spm_bases = {"hrf":[0,0]}

parfile_template = "%s/parfiles/TestRT_r%d_%s_%s.txt"

names = ["simple", "choice", "inst"]

cont01 = ["simple", "T", names, [1,0,0]]
cont02 = ["choice", "T", names, [0,1,0]]
cont03 = ["simple-choice", "T", names, [1,-1,0]]
cont04 = ["choice-simple", "T", names, [-1,1,0]]
