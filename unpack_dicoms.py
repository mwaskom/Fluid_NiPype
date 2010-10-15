"""
DICOM conversion and preprocessing script for Fluid Intelligence project.

Steps:
 - Fetch DICOM files from Sigma
 - Parse DICOM directory to obtain sequence info file
 - Automatically unpack the full acquisitions
 - Create heuristic softlinks to Nifti and MGZ files
 - Run the Freesurfer reconstruction pipeline
 - Preprocesses DWI images
 - Perform spatial normalization to FSL's template with FLIRT and FNIRT
"""
import os
import sys
import shutil
import time
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
parser.add_argument("-fetch", action="store_true",
                    help="run fetch_dicoms to copy from sigma before unpacking")
parser.add_argument("-moco", action="store_true",
                    help="unpack the MoCo BOLD runs")
parser.add_argument("-norun", dest="run", action="store_false",
                    help="don't run the unpacking pipeline")
parser.add_argument("-nolink", dest="link", action="store_false",
                    help="don't create the heuristic links")
parser.add_argument("-norecon", dest="recon", action="store_false",
                    help="Submit a recon-all job to SGE after unpacking")
parser.add_argument("-nodwi", dest="dwi", action="store_false",
                    help="don't run dt_recon to unpack the DWI image")
parser.add_argument("-noreg", dest="reg", action="store_false",
                    help="don't perform FSL registration")
parser.add_argument("-reparse", action="store_true",
                    help="force rerunning of the dicom directory parsing")
parser.add_argument("-reconvert", action="store_true",
                    help="force rerunning of the dicom conversion")
parser.add_argument("-relink", action="store_true",
                    help="overwrite any old heuristic links")
parser.add_argument("-inseries", action="store_true", 
                    help="force running nipype in series")
parser.add_argument("-debug", action="store_true", help="turn on debugging")

# Display help if we got less than one argument
if len(sys.argv) < 2:
    sys.argv.insert(0,"-h")
args = parser.parse_args()

# Hardcoded data directory and template
datadir = "/mindhive/gablab/fluid/Data/"
subject_template = "gf*"

# If type not specified, set to a null string
# Dicom files will be expected to be in $subject/dicom
# This should work, but I haven't actually tested it.
if not args.type:
    args.type = ""

# Don't do any structural processing if we're just 
# unpacking a functional session
if args.type == "func":
    args.recon = False
    args.dwi = False
    args.reg = False

# Initialize the subjectinfo dict
subjectinfo = {}

# Get subjects from command line, or glob based on template
if args.fetch and args.all:
    sys.exit("Cannot use -all when requesting DICOM fetch")
if args.subjects:
    subjects = args.subjects
elif args.all:
    subject_dirs = glob(os.path.join(datadir,subject_template,"dicom",args.type))
    subjects = [d.split("/")[-3] for d in subject_dirs]
else:
    sys.exit("\nMust use either -subjects or -all")
if args.debug:
    print "Subjects: " + " ".join(subjects)

# Fetch the DICOMS
# ----------------
if args.fetch:
    for subject in subjects:
        targdir = os.path.join(datadir, subject, "dicom", args.type)
        os.system("fetch_dicoms -s %s -l -d %s"%(subject, targdir))

# Pipeline functions
# ------------------

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

def fieldmap_type(dcmfile):
    """Determine the type of a fieldmapping dicom"""
    cmd = ["mri_probedicom", "--i", dcmfile, "--t", "8", "8"]
    proc = subprocess.Popen(cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if stdout.strip()[-4] == "M":
        return "mag"
    elif stdout.strip()[-4] == "P":
        return "phase"
    else:
        raise Exception("Could not determine fieldmap type of %s"%dcmfile)

def parse_info_file(dcminfofile, writeflf=True):
    """Get information from the dicom info file about the runs"""
    infopath = dcminfofile.split("/")
    try:
        subject = [d.replace("_sid_","") for d in infopath if d.startswith("_sid_")][0]
    except IndexError:
        subject = [d for d in infopath if d.startswith(subject_template.replace("*",""))][0]
    
    try:
        info = subjectinfo[subject]
        if args.debug:
            print "Using cached info list"
    except KeyError:
        seqfile = np.genfromtxt(dcminfofile, dtype=object)
        info = []
        # Parse by line
        for l in seqfile:
            seqn = l[2]
            if args.debug:
               print "Sequence info:\n\t%s"%(" ".join(l))
            if not args.dontunpack or not seqn in args.dontunpack:
                name = l[12]
                dcmfile = l[1]
                x,y,s,t = tuple([int(dim) for dim in l[6:10]])
                
                if (s==176) and (t==1) and name.startswith("T1_MPRAGE"):
                    info.append(("structural",dcmfile,seqn, name))
                    info.append(("mprage",dcmfile,seqn, name))
                elif name == "ep2d_t1w":
                    if args.type == "full":
                        pfix = "func"
                    else:
                        pfix = args.type
                    info.append(("structural",dcmfile,seqn, "%s-%s"%(pfix,name)))
                elif name.startswith("field_mapping"):
                    fmaptype = fieldmap_type(
                        os.path.join(datadir,subject,"dicom",args.type,dcmfile))
                    if name == "field_mapping":
                        seqname = "func_%s_fm"%fmaptype
                    else:
                        seqname = name.replace("field_mapping_","").lower() + "_%s_fm"%fmaptype
                    info.append(("fieldmaps",dcmfile,seqn,seqname))
                elif (t==70) and name.startswith("DIFFUSION"):
                    info.append(("dwi",dcmfile,seqn,"diffusion"))
                elif (s==176) and (t==1) and name.startswith("gre_mgh_multiecho"):
                    angle = int(name[18:-12])
                    info.append(("flash",dcmfile,seqn,"flash_%02d-rms"%angle))
                elif ("ge_func" in name
                      and not is_moco(os.path.join(datadir,subject,"dicom",args.type,dcmfile))):
                    if (t==198) and name.startswith("NBack"):
                        info.append(("bold",dcmfile,seqn,name))
                    elif (t==248) and name.startswith("MOT_Block"):
                        info.append(("bold",dcmfile,seqn,name))
                    elif (t==302) and name.startswith("MOT_Jitter"):
                        info.append(("bold",dcmfile,seqn,name))
                    elif (t==208) and name.startswith("RT_ge_func"):
                        info.append(("bold",dcmfile,seqn,name))
                    elif (t==244) and name.startswith("IQ_ge_func"):
                        info.append(("bold",dcmfile,seqn,name))
                    elif (t==62) and name.endswith("Resting"):
                        info.append(("bold",dcmfile,seqn,"Resting"))
                elif (args.moco
                      and "ge_func" in name
                      and is_moco(os.path.join(datadir,subject,"dicom",args.type,dcmfile))):
                    name = name.split("_")[0] + "-moco"
                    if (t==198) and name.startswith("NBack"):
                        info.append(("bold",dcmfile,seqn,name))
                    elif (t==248) and name.startswith("MOT_Block"):
                        info.append(("bold",dcmfile,seqn,name))
                    elif (t==302) and name.startswith("MOT_Jitter"):
                        info.append(("bold",dcmfile,seqn,name))
                    elif (t==244) and name.startswith("IQ_ge_func"):
                        info.append(("bold",dcmfile,seqn,name))
                    elif (t==62) and name.endswith("Resting"):
                        info.append(("bold",dcmfile,seqn,"Resting"))
                else:
                    if args.debug:
                        print "***Skipping sequence***"
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


# Unpacking pipeline 
# ------------------

# Subject id iterable
subjsource = pe.Node(interface=util.IdentityInterface(fields=["sid"]),
                     name='subjinfo',
                     iterables = ("sid", subjects))

# Parse a dicom directory to get sequence information
dcminfo = pe.Node(interface=fs.ParseDICOMDir(sortbyrun=True,
                                             summarize=True),
                  name="dicominfo")
dcminfo.overwrite = args.reparse

# DataGrabber for dicom images
dcmsource = pe.MapNode(interface=nio.DataGrabber(infields=["sid","dcmfile"],
                                                 outfields=["dicompath"]),
                       name="dicomsource",
                       iterfield=["dcmfile"])

dcmsource.inputs.base_directory = datadir
dcmsource.inputs.template = "%s/%s/"+args.type+"/%s"
dcmsource.inputs.template_args = dict(dicompath=[["sid","dicom","dcmfile"]])

# DataGrabber for filelist files
flfsource = pe.MapNode(interface=nio.DataGrabber(infields=["sid","seqn"],
                                                 outfields=["flfpath"]),
                       name="flfsource",
                       iterfield=["seqn"])

flfsource.inputs.base_directory = datadir
flfsource.inputs.template = "%s/unpack/flf/%s_"+args.type+"_seq%s"
flfsource.inputs.template_args = dict(flfpath=[["sid","sid","seqn"]])

# Dummy node to facilitate sensible naming
nameimg = pe.MapNode(interface=util.IdentityInterface(fields=["outfile"]),
                     name="nameimage",
                     iterfield=["outfile"])

# Use mri_convert to unpack the dicoms
convert = pe.MapNode(interface=fs.MRIConvert(in_type="siemens_dicom"),
                     name="convertdicoms",
                     iterfield=["in_file", "sdcm_list", "out_type"])
convert.overwrite = args.reconvert

# Node to create a substitutions list
rename = pe.Node(interface=util.Merge(3, axis="hstack"), name="renameniftis")

# DataSink for the nifti/mgz files
datasink = pe.Node(interface=nio.DataSink(),name="datasink")
datasink.inputs.base_directory = datadir
datasink.inputs.parameterization = False

# Workflow definition
unpack = pe.Workflow(name="unpack_%s"%args.type)
unpack.base_dir = "/mindhive/gablab/fluid/Analysis/NiPype/workingdir/unpack"

# Workflow connection
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
    unpack.run(inseries=args.inseries)

    # Log archiving
    logdir = "/mindhive/gablab/fluid/NiPype_Code/log_archive/%s"%datestamp
    if not os.path.isdir(logdir):
        os.mkdir(logdir)
    timestamp = str(datetime.now())[11:16].replace(":","-")
    fulllog = open("%s/unpack_%s.log"%(logdir,timestamp),"w")
    for lf in ["pypeline.log%s"%n for n in [".4",".3",".2",".1",""]]:
        if os.path.isfile(lf):
            fulllog.write(open(lf).read())
            os.remove(lf)
    fulllog.close()
    
    # Write the pipeline graph
    unpack.write_graph(graph2use="flat")


# Softlink heuristically 
# ----------------------

for subj in subjects:
    # This file should be created by the pipeline above
    infofile = os.path.join(datadir, subj, "unpack/%s-dicominfo.txt"%args.type)
    if args.debug:
        print "Parsing %s for linking"%infofile
    if os.path.exists(infofile) and args.link:
        # Get information about the sequences
        info = parse_info_file(infofile, writeflf=False)
        # Hash dict to control names for multiple runs
        boldhash = dict(IQ=1,NBack=1,RT=1,MOT_Jitter=1,MOT_Block=1)
        for seq in info:
            if args.debug:
                print "Sequence info: %s"%str(seq)
            # This is the recon-all mprage souce
            if seq[0]=="mprage":
                src = get_img_name(subj, seq[2], "mgz")
                trg = os.path.join(datadir, subj, "mri/orig/001.mgz")
            # Other images get linked here
            else:
                type = seq[0]
                dcmfile = seq[1]
                seqn = seq[2]
                name = seq[3]
                trgdir = os.path.join(datadir, subj, type)
                # Nifti mprage, low-res T1, and FLASH images
                if type == "structural":
                    if name.endswith("ep2d_t1w"):
                        src = get_img_name(subj, seqn, "nii.gz")
                        trgfile = name + ".nii.gz"
                    elif name.startswith("T1_MPRAGE"):
                        src = get_img_name(subj, seqn, "nii.gz")
                        trgfile = "mprage.nii.gz"
                # We unpack the FLASH rms image here, and get
                # the individual echos later
                elif type == "flash":
                    src = get_img_name(subj, seqn, "mgz")
                    trgfile = "%s.mgz"%name
                # Functional runs
                elif type == "bold":
                    src = get_img_name(subj, seqn, "nii.gz")
                    if is_moco(dcmfile):
                        moco = "moco"
                    else:
                        moco = ""
                    if name == "Resting":
                        trgfile = "Resting%s.nii.gz"%moco
                    else:
                        if name.startswith("MOT"):
                            par = "MOT_" + name.split("_")[1]
                        else:
                            par = name.split("_")[0]
                        nrun = boldhash[par]
                        boldhash[par] += 1
                        trgfile = "%s_run%d%s.nii.gz"%(par,nrun,moco)
                # Diffusion.  Not calling it DTI, to please Satra
                elif type == "dwi":
                    src = get_img_name(subj, seqn, "nii.gz")
                    trgfile = "%s.nii.gz"%name
                # Fieldmaps.  Both magnitude and phase series get unpacked
                elif type == "fieldmaps":
                    src = get_img_name(subj, seqn, "nii.gz")
                    trgfile = "%s.nii.gz"%name
                else:
                    pass
                trg = os.path.join(trgdir, trgfile)
            if os.path.exists(src):
                if os.path.lexists(trg) and args.relink:
                    os.remove(trg)
                if not os.path.lexists(trg):
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
    elif args.link:
        print "Info file %s not found"%infofile

# Structural Preprocessing
# ------------------------

for subj in subjects:
    sgescript = []

    # Run recon-all 
    # ---------------
    if args.recon and os.path.exists(os.path.join(datadir, subj, "mri/orig/001.mgz")):
        # Make sure a recon hasn't been started for this subject
        if not os.path.isfile(os.path.join(datadir, subj, "scripts/recon-all-status.log")):
            # Recon-all command line
            sgescript.append("recon-all -s %s -all"%subj)
            print "Adding %s recon job to SGE script"%subj
        else:
            print "Recon submission requested for %s, but recon status log exists"%subj
    elif args.recon:
        print "Recon submission requested for %s, but recon source image does not exist"%subj
    
    # Unpack the DWI image with dt_recon
    # ----------------------------------
    infofile = os.path.join(datadir, subj, "unpack/%s-dicominfo.txt"%args.type)
    if args.dwi and os.path.exists(infofile):
        # Figure out of there's a diffusion acquisision in our dicoms
        diffinfo = [i for i in parse_info_file(infofile, writeflf=False) if i[0]=="dwi"][0]
        if diffinfo:
            if args.debug:
                print diffinfo
            srcfile = os.path.join(datadir, subj, "dicom", args.type, diffinfo[1])
            # Make sure srcfile exists
            if os.path.exists(srcfile):
                trgdir = os.path.join(datadir, subj, "dwi")
                if not os.path.exists(trgdir):
                    os.makedirs(trgdir)
                # Check for this process by looking for the log file
                if not os.path.exists(os.path.join(trgdir, "dt_recon.log")):
                    # DT_recon command line
                    sgescript.append("dt_recon --i %s --s %s --o %s"%(srcfile, subj, trgdir))
                    print "Adding %s dt_recon job to SGE script"%subj
                else:
                    print "dt_recon log found for %s; skipping dwi unpacking"%subj
            else:
                "DWI source DICOM not found for %s, skipping dwi unpacking"%subj

    # Perform spatial normalization
    # -----------------------------
    if args.reg and (os.path.exists(os.path.join(datadir, subj, "mri", "brainmask.mgz")) and
                     os.path.exists(os.path.join(datadir, subj, "mri", "T1.mgz"))):
        # Fluid_register command line
        sgescript.append(
            "python /mindhive/gablab/fluid/NiPype_Code/fluid_normalize.py %s"%subj)
        print "Adding %s normalization to SGE script"%subj
    elif args.reg:
        print "Normalization requested for %s, but source images not found"%subj
        
    # Actually submit to Sun Grid Engine
    # ----------------------------------
    if sgescript:
        sgedir = os.path.join(os.getenv("HOME"), "sge")
        if not os.path.exists(sgedir):
            os.mkdir(sgedir)
        # Write a script to give to qsub
        scriptfile = os.path.join(sgedir, "%s_recon_%s.sh"%(subj, time.time()))
        fid = open(scriptfile,"w")
        fid.write("\n".join(["#! /bin/bash",
                             "#$ -cwd",
                             "#$ -N %s_recon"%subj,
                             "#$ -r y",
                             "#$ -S /bin/bash",
                             ""]))
        fid.write("\n".join(sgescript))
        fid.close()
        
        # Write the qsub command line
        qsub = ["cd",sgedir,";","qsub","-q","long",scriptfile]
        cmd = " ".join(qsub)
        
        # Submit the job
        print "Submitting job %s_recon to Sun Grid Engine"%subj
        proc = subprocess.Popen(cmd,
                                stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE,
                                env=os.environ,
                                shell=True,
                                cwd=sgedir)
        stdout, stderr = proc.communicate()
        
