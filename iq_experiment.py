"""
    Nipype experiment module for the NBack paradigm.
"""
template_args = dict(func=[["subject_id", "IQ*_run1"]],
                    target=[["subject_id", "ep2d_t1w"]],
                    struct=[["subject_id", "mprage"]])

source_template = "%s/nii/%s.nii.gz"

hpcutoff = 128
TR = 2.
units = "secs"

nruns = 1

exclude_subjects = ["SMARTER_SP%02d"%i for i in [2,4] + range(7,21)] 

fsl_bases = {"dgamma":{"derivs":False}}
spm_bases = {"hrf":[0,0]}

parfile_template = "%s/parfiles/IQ_scan1_r%d_%s_%s.txt"

names = ["easy", "hard"]

cont01 = ["easy", "T", names, [1,0]]
cont02 = ["hard", "T", names, [0,1]]
cont03 = ["easy-hard", "T", names, [1,-1]]
cont04 = ["hard-easy", "T", names, [-1,1]]
