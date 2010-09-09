"""
DICOM Unpacking script for the Fluid Intelligence project.
Will unpack DICOM directory, create heuristic softlinks, and optionally
submit recon-all jobs to the Sun Grid Engine
"""
import os
import sys
import shutil
from datetime import datetime
import subprocess
import argparse
from glob import glob

import numpy as np
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
import nipype.interfaces.freesurfer as fs
import nipype.pipeline.engine as pe

# Command line arguments
parser = argparse.ArgumentParser(usage="unpack_dicoms.py [options]")
parser.add_argument("-subjects", nargs="*", metavar="subjid",
                    help="list of subject ids to unpack")
parser.add_argument("-type", metavar="scantype", help="func, struct, or full")
parser.add_argument("-all", action="store_true", 
                    help="run all subjects with dicom dir")
parser.add_argument("-dontunpack", nargs="*", metavar="num",
                    help="sequence number of run(s) to skip")
parser.add_argument("-moco", action="store_true",
                    help="unpack the MoCo BOLD runs")
parser.add_argument("-recon", action="store_true",
                    help="Submit a recon-all job to SGE after unpacking")
parser.add_argument("-norun", dest="run", action="store_false",
                    help="don't run the unpacking pipeline")
parser.add_argument("-nolink", dest="link", action="store_false",
                    help="don't create the heuristic links")
parser.add_argument("-reparse", action="store_true",
                    help="force rerunning of the dicom directory parsing")
parser.add_argument("-reconvert", action="store_true",
                    help="force rerunning of the dicom conversion")
parser.add_argument("-relink", action="store_true",
                    help="overwrite any old heuristic links")
parser.add_argument("-inseries", action="store_true", 
                    help="force running nipype in series")
parser.add_argument("-debug", action="store_true", help="turn on debugging")
args = parser.parse_args()

# Hardcoded data directory and template
datadir = "/mindhive/gablab/fluid/Data/"
subject_template = "gf*"

# Initialize the subjectinfo dict
subjectinfo = {}

# Get subjects from command line, or glob based on template
if args.subjects:
    subjects = args.subjects
elif args.all:
    subject_dirs = glob(os.path.join(datadir, subject_template, "dicom", args.type))
    subjects = [d.split("/")[-2] for d in subject_dirs]
else:
    sys.exit("\nMust use either -subjects or -all")

#====================#
# Pipeline functions #
#====================#

def get_dicom_dir(subject):
    """Return the path to a dicom directory"""
    return glob(os.path.join(datadir,subject,"dicom",args.type))[0]

def is_moco(dcmfile):
    """Determine if a run has on-line motion correction"""
    cmd = ['mri_probedicom', '--i', dcmfile, '--t', '8', '103e']
    proc  = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    return stdout.strip().startswith('MoCoSeries')

def parse_info_file(dcminfofile, writeflf=True):
    """Get information from the dicom info file about the runs"""
    infopath = dcminfofile.split("/")
    try:
        subject = [d.replace("_sid_","") for d in infopath if d.startswith("_sid_")][0]
    except IndexError:
        subject = [d for d in infopath if d.startswith(subject_template.replace("*",""))][0]
    
    try:
        info = subjectinfo[subject]
    except KeyError:
        seqfile = np.genfromtxt(dcminfofile, dtype=object)
        info = []
        # Parse by line
        for l in seqfile:
            seqn = l[2]
            if not args.dontunpack or not seqn in args.dontunpack:
                name = l[12]
                dcmfile = l[1]
                x,y,s,t = tuple([int(dim) for dim in l[6:10]])
                
                if (s==176) and (t==1) and name.startswith("T1_MPRAGE"):
                    info.append(("structural",dcmfile,seqn, name))
                    info.append(("mprage",dcmfile,seqn, name))
                elif (t==8) and name.startswith("gre_mgh_multiecho"):
                    angle = name[18:-12]
                    info.append(("structural",dcmfile,seqn, "flash_%02d"%int(angle)))
                elif name == "ep2d_t1w":
                    info.append(("structural",dcmfile,seqn, name))
                elif not is_moco(os.path.join(datadir,subject,"dicom",args.type,dcmfile)):
                    if (t==198) and name.startswith("NBack"):
                        info.append(("bold",dcmfile,seqn, name))
                    elif (t==208) and name.startswith("RT_ge_func"):
                        info.append(("bold",dcmfile,seqn, name))
                    elif (t==244) and name.startswith("IQ_ge_func"):
                        info.append(("bold",dcmfile,seqn, name))
                    elif (t==62) and name.endswith("Resting"):
                        info.append(("bold",dcmfile,seqn,"Resting"))
                elif (args.moco
                      and is_moco(os.path.join(datadir,subject,"dicom",args.type,dcmfile))):
                    name = name.split("_")[0] + "-moco"
                    if (t==198) and name.startswith("NBack"):
                        info.append(("bold",dcmfile,seqn, name))
                    elif (t==208) and name.startswith("RT_ge_func"):
                        info.append(("bold",dcmfile,seqn, name))
                    elif (t==244) and name.startswith("IQ_ge_func"):
                        info.append(("bold",dcmfile,seqn, name))
                    elif (t==62) and name.endswith("Resting"):
                        info.append(("bold",dcmfile,seqn,"Resting"))
                else:
                    pass
        subjectinfo[subject] = info

    # Write files containing lists of the dicoms in each series to speed up mri_convert
    if writeflf:
        for seqinfo in info:
            flfdir = os.path.join(datadir, subject, "unpack", "flf")
            if not os.path.isdir(flfdir):
                os.makedirs(flfdir)
            alldcms = glob(os.path.join(
                datadir,subject,"dicom",args.type,seqinfo[1][:-5]+"*"))
            flffile = open(os.path.join(
                flfdir, subject + "_%s"%args.type + "_seq" + seqinfo[2]),"w")
            flffile.write(" ".join(alldcms))
            flffile.close()

    return info

def get_out_filename(dcminfofile):
    """Name the raw nifti files based on the seqence number"""
    infolist = parse_info_file(dcminfofile)
    return ["%s_seq_%02d"%(args.type,int(t[2])) for t in infolist]

def get_out_ftype(dcminfofile):
    """Write files to nii.gz unless it's the mprage or flash"""
    infolist = parse_info_file(dcminfofile)
    typelist = []
    if args.debug:
        print "Running get_out_ftype function"
    for info in infolist:
        if args.debug:
            print "Sequence info: %s"%str(info)
        if info[0] == "mprage" or info[3].startswith("flash"):
            typelist.append(("mgz", ".mgz"))
        else:
            typelist.append(("niigz", ".nii.gz"))
    return typelist
    
def get_img_name(subj, seq, ftype):
    """Get the source file name"""    
    return os.path.join(datadir, subj, "nifti/%s_seq_%02d.%s"%(args.type,int(seq),ftype))


#====================#
# Unpacking pipeline #
#====================#

subjsource = pe.Node(interface=util.IdentityInterface(fields=["sid"]),
                     name='subjinfo',
                     iterables = ("sid", subjects))

dcminfo = pe.Node(interface=fs.ParseDICOMDir(sortbyrun=True,
                                             summarize=True),
                  name="dicominfo")
dcminfo.overwrite = args.reparse

dcmsource = pe.MapNode(interface=nio.DataGrabber(infields=["sid","dcmfile"],
                                                 outfields=["dicompath"]),
                       name="dicomsource",
                       iterfield=["dcmfile"])

dcmsource.inputs.base_directory = datadir
dcmsource.inputs.template = "%s/%s/"+args.type+"/%s"
dcmsource.inputs.template_args = dict(dicompath=[["sid","dicom","dcmfile"]])

flfsource = pe.MapNode(interface=nio.DataGrabber(infields=["sid","seqn"],
                                                 outfields=["flfpath"]),
                       name="flfsource",
                       iterfield=["seqn"])

flfsource.inputs.base_directory = datadir
flfsource.inputs.template = "%s/unpack/flf/%s_"+args.type+"_seq%s"
flfsource.inputs.template_args = dict(flfpath=[["sid","sid","seqn"]])

nameimg = pe.MapNode(interface=util.IdentityInterface(fields=["outfile"]),
                     name="nameimage",
                     iterfield=["outfile"])

convert = pe.MapNode(interface=fs.MRIConvert(in_type="siemens_dicom"),
                     name="convertdicoms",
                     iterfield=["in_file", "sdcm_list", "out_type"])
convert.overwrite = args.reconvert

rename = pe.Node(interface=util.Merge(3, axis="hstack"), name="renameniftis")

datasink = pe.Node(interface=nio.DataSink(),name="datasink")
datasink.inputs.base_directory = datadir
datasink.inputs.parameterization = False

unpack = pe.Workflow(name="unpackdicoms")
unpack.base_dir = "/mindhive/gablab/fluid/Analysis/NiPype/workingdir/unpack"

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
     (dcminfo, convert,
        [(("dicom_info_file", lambda x: [t[0] for t in get_out_ftype(x)]), "out_type")]),
     (dcmsource, convert,
        [("dicompath", "in_file")]),
     (flfsource, convert,
        [("flfpath", "sdcm_list")]),
     (convert, rename,
        [(("out_file", lambda x: [os.path.split(f)[1] for f in x]), "in1")]),
     (nameimg, rename,
        [("outfile", "in2")]),
     (dcminfo, rename,
        [(("dicom_info_file", get_out_ftype), "in3")]),
     (subjsource, datasink,
        [("sid", "container")]),
     (rename, datasink,
        [(("out", lambda x: 
            [("dicominfo.txt","%s-dicominfo.txt"%args.type)] + 
                [(l[0],l[1]+l[2][1]) for l in x]), "substitutions")]),
     (convert, datasink,
        [("out_file", "nifti.@niifiles")]),
     (dcminfo, datasink,
        [("dicom_info_file", "unpack.@info")]),
    ])

# File crashdumps by date
datestamp = str(datetime.now())[:10]
codepath = os.path.split(os.path.abspath(__file__))[0]
crashdir = os.path.abspath("%s/crashdumps/%s" % (codepath, datestamp))
if not os.path.isdir(crashdir):    
    os.makedirs(crashdir)
unpack.config = dict(crashdump_dir=crashdir) 

# Run the pipeline
if args.run:
    timestamp = str(datetime.now())[:-10].replace("-","").replace(":","").replace(" ","-")
    unpack.run(inseries=args.inseries)
    fulllog = open("/mindhive/gablab/fluid/NiPype_Code/log_archive/unpack-%s.log"%timestamp,"w")
    for lf in ["pypeline.log%s"%n for n in [".4",".3",".2",".1",""]]:
        if os.path.isfile(lf):
            fulllog.write(open(lf).read())
            os.remove(lf)
    fulllog.close()
    unpack.write_graph(graph2use="flat")

#========================#
# Softlink heuristically #
#========================#

for subj in subjects:
    infofile = os.path.join(datadir, subj, "unpack/%s-dicominfo.txt"%args.type)
    if args.debug:
        print "Parsing %s for linking"%infofile
    if os.path.isfile(infofile) and args.link:
        info = parse_info_file(infofile, writeflf=False)
        boldhash = dict(IQ=1,NBack=1,RT=1)
        for seq in info:
            if args.debug:
                print "Sequence info: %s"%str(seq)
            if seq[0]=="mprage":
                src = get_img_name(subj, seq[2], "mgz")
                trg = os.path.join(datadir, subj, "mri/orig/001.mgz")
            else:
                type = seq[0]
                name = seq[3]
                trgdir = os.path.join(datadir, subj, type)
                if type == "structural":
                    if name == "ep2d_t1w":
                        src = get_img_name(subj, seq[2], "nii.gz")
                        trgfile = name + ".nii.gz"
                    elif name.startswith("T1_MPRAGE"):
                        src = get_img_name(subj, seq[2], "nii.gz")
                        trgfile = "mprage.nii.gz"
                    elif name.startswith("flash"):
                        src = get_img_name(subj, seq[2], "mgz")
                        trgfile = "%s.mgz"%name
                elif type == "bold":
                    src = get_img_name(subj, seq[2], "nii.gz")
                    if name == "Resting":
                        trgfile = "Resting.nii.gz"
                    else:
                        par = name.split("_")[0]
                        nrun = boldhash[par]
                        boldhash[par] += 1
                else:
                    pass
                trg = os.path.join(trgdir, trgfile)
            if os.path.isfile(src):
                if os.path.isfile(trg) and args.relink:
                    os.remove(trg)
                if not os.path.isfile(trg):
                    trgdir = os.path.split(trg)[0]
                    if not os.path.isdir(trgdir):
                        os.makedirs(trgdir)
                    if args.debug:
                        print "Linking %s to %s"%(src, trg)
                    os.symlink(src, trg)
                else:
                    print "Target file %s exists; use -relink to overwrite"%trg
            else:
                print "Source file %s not found"%src
    elif not args.link:
        print "Info file %s not found"%infofile

#=================#
# Start recon-all #
#=================#
if args.recon:
    for subj in subjects:
        # Make sure the source image exists
        if os.path.isfile(os.path.join(datadir, subj, "mri/orig/001.mgz")):
            # Make sure a recon hasn't been started for this subject
            if not os.path.isfile(os.path.join(datadir, subj, "scripts/recon-all-status.log")):
                reconcmd = "'recon-all -s %s -all'"%subj
                sgecmd = "ezsub.py -c %s -n %s_recon -q long.q"%(reconcmd, subj)
                print "Submitting %s recon job to SGE"%subj
                os.system(sgecmd)
            else:
                print "Recon submission requested for %s, but recon status log exists"%subj
        else:
            print "Recon submission requested for %s, but recon source image does not exist"%subj
