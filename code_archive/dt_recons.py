"""Preprocess dwi data for the subjects we've already unpacked"""
import os
import sys
import numpy as np
from glob import glob

data_dir = "/mindhive/gablab/fluid/Data"

diffniftis = glob(os.path.join(data_dir, "gf??/dwi/diffusion.nii.gz"))

subjects = [p.split("/")[-3] for p in diffniftis]

def get_dwi_dicom(subj):
    print subj
    for type in ["struct", "full"]:
        dfile = os.path.join(data_dir,subj,"unpack","%s-dicominfo.txt"%type)
        if os.path.exists(dfile):
            seqinfo = np.genfromtxt(dfile,dtype=object)
            scantype = type
    dicoms = []
    for l in seqinfo:
        name = l[-1]
        dcmfile = l[1]
        t = l[9]
        if name.startswith("DIFF") and (t=="70"):
            dicoms.append(dcmfile)
    if len(dicoms) != 1:
        print "Subject %s did not match 1 DWI acquisition"%subj
    else:
        return (dicoms[0], scantype)

subjects.remove("gf15")
for subj in subjects:
    dcmfile, scantype = get_dwi_dicom(subj)
    srcfile = os.path.join(data_dir, subj, "dicom", scantype, dcmfile)
    trgdir = os.path.join(data_dir, subj, "dwi")
    cmd = "dt_recon --i %s --s %s --o %s"%(srcfile, subj, trgdir)
    os.system(cmd)
