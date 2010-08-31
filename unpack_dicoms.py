import os
import sys
import shutil
import subprocess
import argparse
from glob import glob

import numpy as np
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
import nipype.interfaces.freesurfer as fs
import nipype.pipeline.engine as pe

parser = argparse.ArgumentParser()
parser.add_argument("-subjects", nargs="*", help="list of subject ids")
parser.add_argument("-all", action="store_true", 
                    help="run all subjects with dicom dir")
args = parser.parse_args()

datadir = "/mindhive/gablab/fluid/Data/"

if args.subjects:
    subjects = args.subjects
elif args.all:
    subject_template = "gf*"
    subject_dirs = glob(os.path.join(datadir, subject_template, "dicom"))
    subjects = [d.split("/")[-2] for d in subject_dirs]
else:
    print "\nMust use either '-subjects' or '-all'"
    sys.exit(0)


def get_dicom_dir(subject):
    return glob(os.path.join(datadir,subject,"dicom"))[0]

def is_moco(dcmfile):
    cmd = ['mri_probedicom', '--i', dcmfile, '--t', '8', '103e']
    proc  = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    return stdout.strip().startswith('MoCoSeries')

def parse_info_file(dcminfofile):

    seqfile = np.genfromtxt(dcminfofile, dtype=object)
    info = []
    for l in seqfile:
        seqn = l[2]
        name = l[12]
        dcmfile = l[1]
        x,y,s,t = tuple([int(dim) for dim in l[6:10]])
        
        if (s==176) and (t==1) and name.startswith("T1_MPRAGE"):
            info.append(("struct",dcmfile,seqn))
            info.append(("mprage",dcmfile,seqn))
        elif (t==8) and name.startswith("gre_mgh_multiecho"):
            info.append(("struct",dcmfile,seqn))
        elif name == "ep2d_t1w":
            info.append(("struct",dcmfile,seqn))
        elif not is_moco(dcmfile):
            if (t==198) and name.startswith("NBack"):
                info.append(("bold",dcmfile,seqn))
            elif (t==208) and name.startswith("RT_ge_func"):
                info.append(("bold",dcmfile,seqn))
            elif (t==244) and name.startswith("IQ_ge_func"):
                info.append(("bold",dcmfile,seqn))
        else:
            pass

    subject = [d.replace("_sid_","") for d in dcminfofile.split("/") if d.startswith("_sid_")][0]
    for seqinfo in info:
        flfdir = os.path.join(datadir, subject, "unpack", "flf")
        if not os.path.isdir(flfdir):
            os.makedirs(flfdir)
        alldcms = glob(os.path.join(datadir,subject,"dicom",seqinfo[1][:-5]+"*"))
        flffile = open(os.path.join(flfdir, subject + "_seq" + seqinfo[2]),"w")
        flffile.write(" ".join(alldcms))
        flffile.close()

    return info

def get_out_filename(dcminfofile):
    infolist = parse_info_file(dcminfofile)
    namelist = []
    for info in infolist:
        if info[0] == "mprage":
            namelist.append("orig.mgz")
        else:
            namelist.append("seq_%s.nii.gz"%info[2])
    return namelist


subjsource = pe.Node(interface=util.IdentityInterface(fields=["sid"]),
                     name='subjinfo',
                     iterables = ("sid", subjects))

dcminfo = pe.Node(interface=fs.ParseDICOMDir(sortbyrun=True,
                                             summarize=True),
                  name="dicominfo")

dcmsource = pe.MapNode(interface=nio.DataGrabber(infields=["sid","dcmfile"],
                                                 outfields=["dicompath"]),
                       name="dicomsource",
                       iterfield=["dcmfile"])

dcmsource.inputs.base_directory = datadir
dcmsource.inputs.template = "%s/%s/%s"
dcmsource.inputs.template_args = dict(dicompath=[["sid","dicom","dcmfile"]])

flfsource = pe.MapNode(interface=nio.DataGrabber(infields=["sid","seqn"],
                                                 outfields=["flfpath"]),
                       name="flfsource",
                       iterfield=["seqn"])

flfsource.inputs.base_directory = datadir
flfsource.inputs.template = "%s/unpack/flf/%s_seq%s"
flfsource.inputs.template_args = dict(flfpath=[["sid","sid","seqn"]])

nameimg = pe.MapNode(interface=util.IdentityInterface(fields=["outfile"]),
                     name="nameimage",
                     iterfield=["outfile"])

convert = pe.MapNode(interface=fs.MRIConvert(in_type="siemens_dicom"),
                     name="convertdicoms",
                     iterfield=["in_file", "out_file", "sdcm_list"])

datasink = pe.Node(interface=nio.DataSink(),name="datasink")
datasink.inputs.base_directory = datadir
datasink.inputs.parameterization = False

unpack = pe.Workflow(name="unpackdicoms")
unpack.base_dir = os.path.join(datadir,"unpack_workingdir")

unpack.connect(
    [(subjsource, dcminfo,
        [(("sid", get_dicom_dir), "dicom_dir")]),
     (subjsource, dcmsource,
        [("sid", "sid")]),
     (dcminfo, dcmsource,
        [(("dicom_info_file", lambda x: [t[1] for t in parse_info_file(x)]),"dcmfile")]),
     (subjsource, flfsource,
        [("sid", "sid")]),
     (dcminfo, flfsource,
        [(("dicom_info_file", lambda x: [t[2] for t in parse_info_file(x)]),"seqn")]),
     (dcminfo, nameimg,
        [(("dicom_info_file", get_out_filename), "outfile")]),
     (dcmsource, convert,
        [("dicompath", "in_file")]),
     (flfsource, convert,
        [("flfpath", "sdcm_list")]),
     (nameimg, convert,
        [("outfile", "out_file")]),
     (subjsource, datasink,
        [("sid", "container")]),
     (convert, datasink,
        [("out_file", "nii.@niftigz")]),
     (dcminfo, datasink,
        [("dicom_info_file", "unpack.@info")]),
    ])

if __name__ == "__main__":
    unpack.run()
    unpack.write_graph(graph2use="flat")
    
