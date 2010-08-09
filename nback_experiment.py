# vi: set ft=python sts=4 ts=4 sw=4 et:
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
"""
    Nipype experiment module for the NBack paradigm.
"""

datainfo = dict(func=[["subject_id", "nii", "nback_r2"]], #,"nback_r2"]]],
            struct=[["subject_id", "nii", "mprage_mc"]])
hpcutoff = 128
TR = 2.

nruns = 1

bases = {"dgamma":{"derivs":False}}

names = ["0back","1back","2back","4back","inst"]
partemp = "%s/parfiles/par_%s_d1_r%d_%s.txt"

cont01 = ["0Back", "T", names, [1,0,0,0,0]]
cont02 = ["1Back", "T", names, [0,1,0,0,0]]
cont03 = ["2Back", "T", names, [0,0,1,0,0]]
cont04 = ["4Back", "T", names, [0,0,0,1,0]]
cont05 = ["1B-0B", "T", names, [-1,1,0,0,0]]
cont06 = ["2B-0B", "T", names, [-1,0,1,0,0]]
cont07 = ["2B-1B", "T", names, [0,-1,1,0,0]]
cont08 = ["4B-0B", "T", names, [-1,0,0,1,0]]
cont09 = ["4B-1B", "T", names, [0,-1,0,1,0]]
cont10 = ["4B-2B", "T", names, [0,0,-1,1,0]]
cont11 = ["421B-0B", "T", names, [-1,1./3,1./3,1./3,0]]
