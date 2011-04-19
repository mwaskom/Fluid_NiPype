#! /usr/bin/env python
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
import time
import shutil
from os.path import join as pjoin
import subprocess
import argparse
from glob import glob
import numpy as np

import nipype.interfaces.utility as util
import nipype.interfaces.io as io
import nipype.interfaces.freesurfer as fs
import nipype.pipeline.engine as pe

from fluid_utility import archive_crashdumps

def main(arglist):

    data_dir = os.environ["DATA"]

    # Parse the commandline
    args = parse_cmdline(arglist)

    # Set all if no processing stages were specified
    doall = False
    if args.procs is None:
        doall = True

    # We'll refer to subject list a lot, so get a direct ref
    s_list = args.subjects

    # Fetch the dicoms from the sigma server
    if doall or "fetch" in args.procs:
        types = fetch_dicoms(data_dir, s_list, args.type)
    elif args.type is None:
        sys.exit("Scan type must be specified if skipping dicom fetch")
    else:
        types = [args.type for s in s_list]

    # Use Nipype to unpack the DICOMS
    if doall or "unpack" in args.procs:
        # Get the unpacking workflow
        unpacker = get_unpack_workflow(s_list)

        # Set the base directory
        # XXX running in /tmp doesn't work with IPython but let's make it
        # so that this gets removed
        unpacker.base_dir = pjoin(data_dir, "unpacking")

        # Feed some information into the workflow
        unpacker.inputs.flowinfo.data_dir = data_dir
        unpacker.inputs.flowinfo.dont_unpack = args.dontunpack
        unpacker.inputs.flowinfo.subj_dict = dict(zip(s_list, types))

        # Set the datasink output directory
        unpacker.inputs.datasink.base_directory = data_dir
        
        # Possibly set the working directory
        if args.workingdir is not None:
            unpacker.base_dir = args.workingdir

        # Archive crashdumps
        archive_crashdumps(unpacker)
        
        # Run the workflow
        if args.ipython:
            plugin = "IPython"
        else:
            plugin = "Linear"
        unpacker.run(plugin=plugin)
    
    if doall or "cache" in args.procs:
        # Cache the runinfo to a .npy file
        cache_run_info(data_dir, s_list, types, args.dontunpack)

    if doall or "link" in args.procs:
        # Create heuristic symlinks to the converted files
        make_symlinks(data_dir, s_list, types, args.dontunpack)

    # Figure out that anatomical preprocessing that needs to be done
    anat_stages = ["recon", "dwi", "norm"]

    if doall or "anat" in args.procs:
        anat_proc = anat_stages
    else:
        anat_proc = []
        for proc in args.procs:
            if proc in anat_stages:
                anat_proc.append(proc)
   
    # And do it
    # (This ends up writing a script and submitting it to the
    # Sun Grid Engine via qsub, but that all happens below)
    if anat_proc:
        preprocess_anatomicals(data_dir, s_list, types, anat_proc)

    # Tar and Gzip the DICOMS
    if doall or "zip" in args.procs:
        compress_dicoms(data_dir, s_list, types)


def fetch_dicoms(data_dir, subjects, type=None):
    """Get the source images and determine scan type from number of dicom files."""
    
    types = []
    
    # Without a given type, sync to temporary directory 
    if type is None:
        syncdir = "sync"
    else:
        syncdir = type

    # Make sure subjects is a list to avoid bad things
    if not isinstance(subjects, list):
        subjects = [subjects]

    # Iterate through the list of subjects, grabbing their dicoms
    for subj in subjects:
        dicomdir = pjoin(data_dir, subj, "dicom")
        targdir = pjoin(dicomdir, syncdir)

        # Call external python script to copy over dicom files from sigma
        fetch_cmd = "fetch_dicoms -s %s -l -q -d %s"%(subj, targdir)
        os.system(fetch_cmd)

        # Figure out scan type based on number of dicom files
        if type is None:
            nfiles = len(glob(pjoin(targdir, "*.dcm")))
            if nfiles < 5000:
                stype = "func"
            elif nfiles < 8500:
                stype = "struct"
            else:
                stype = "full"
            os.rename(targdir, pjoin(dicomdir, stype))
        else:
            stype = type
        types.append(stype)

    # Return a list of scan type for each subject
    return types

def get_unpack_workflow(subject_list=[]):
    """Define a nipype workflow to unpack the dicoms."""

    # Create the workflow
    unpacker = pe.Workflow(name="unpack")

    # Set up some informational nodes
    subjsource = pe.Node(util.IdentityInterface(fields=["subject_id"]),
                         iterables=("subject_id", subject_list),
                         name="subjsource")
    

    flowinfo = pe.Node(util.Function(input_names=["data_dir",
                                                  "dont_unpack",
                                                  "subj_dict",
                                                  "subject_id"],
                                     output_names=["dicom_dir",
                                                   "dont_unpack",
                                                   "scan_type"],
                                     function=flow_info),
                         name="flowinfo")
    
    # Parse the dicom directory to get sequence information
    parsedcms = pe.Node(interface=fs.ParseDICOMDir(sortbyrun=True,
                                                   summarize=True),
                        name="parsedcms")

    # Mine the dicom directory summary file for information about runs
    runinfo = pe.Node(interface=util.Function(input_names=["info_file",
                                                           "dicom_dir",
                                                           "dont_unpack"],
                                              output_names=["files",
                                                            "groups",
                                                            "seqs",
                                                            "names",
                                                            "ftypes",
                                                            "args"],
                                              function=run_info),
                      name="runinfo")

    # Get a list of dicom files we want to convert
    dcmpath = pe.Node(util.Function(input_names=["dicom_dir",
                                                 "dicom_files",
                                                 "args"],
                                    output_names=["out_files"],
                                    function=get_dicom_path),
                       name="dcmpath")

    # Write text files identifying the dicoms in a series for faster processing
    writeflfs = pe.Node(util.Function(input_names=["dicom_dir",
                                                   "dicom_files",
                                                   "seq_numbers",
                                                   "groups",
                                                   "args"],
                                     output_names=["fl_files"],
                                     function=write_file_lists),
                        name="writeflfs")

    # Convert the files from DCM
    convert = pe.MapNode(fs.MRIConvert(in_type="siemens_dicom"),
                         iterfield=["in_file", "sdcm_list", "out_type", "args"],
                         name="convert")
                                       
    # Rename the image files according to sequence
    rename = pe.Node(util.Function(input_names=["scan_type",
                                                "converted_images",
                                                "seq_numbers",
                                                "ftypes",
                                                "args"],
                                   output_names=["out_files"],
                                   function=rename_images),
                     name="rename")

    # Rename the info file according to scan type
    renameinfo = pe.Node(util.Function(input_names=["scan_type",
                                                    "info_file"],
                                       output_names=["info_file"],
                                       function=rename_info_file),
                         name="renameinfo") 

    # Sink the files in the data directory
    datasink = pe.Node(io.DataSink(parameterization=False),
                       name="datasink")

    unpacker.connect([
        (subjsource,  flowinfo,     [("subject_id", "subject_id")]),
        (flowinfo,    parsedcms,    [("dicom_dir", "dicom_dir")]),
        (flowinfo,    runinfo,      [("dont_unpack", "dont_unpack"),
                                     ("dicom_dir", "dicom_dir")]),
        (parsedcms,   runinfo,      [("dicom_info_file", "info_file")]),
        (flowinfo,    dcmpath,      [("dicom_dir", "dicom_dir")]),
        (runinfo,     dcmpath,      [("files", "dicom_files"),
                                     ("args","args")]),
        (flowinfo,    writeflfs,    [("dicom_dir", "dicom_dir")]),
        (runinfo,     writeflfs,    [("files", "dicom_files"),
                                     ("seqs", "seq_numbers"),
                                     ("groups", "groups"),
                                     ("args", "args")]),
        (dcmpath,     convert,      [("out_files", "in_file")]),
        (writeflfs,   convert,      [("fl_files", "sdcm_list")]),
        (runinfo,     convert,      [("ftypes", "out_type"),
                                     ("args","args")]),
        (flowinfo,    rename,       [("scan_type", "scan_type")]),
        (convert,     rename,       [("out_file", "converted_images")]),
        (runinfo,     rename,       [("seqs", "seq_numbers"),
                                     ("ftypes", "ftypes"),
                                     ("args", "args")]),
        (parsedcms,   renameinfo,   [("dicom_info_file", "info_file")]),
        (flowinfo,    renameinfo,   [("scan_type", "scan_type")]),
        (rename,      datasink,     [("out_files", "nifti.@images")]),
        (renameinfo,  datasink,     [("info_file", "unpack.@info")]),
        (subjsource,  datasink,     [("subject_id", "container")]),
    ])

    return unpacker

def flow_info(data_dir, dont_unpack, subj_dict, subject_id):
    """Little function to provide some basic info to the workflow.""" 
    from os.path import join as pjoin

    scan_type = subj_dict[subject_id]
    dicom_dir = pjoin(data_dir, subject_id, "dicom", scan_type)
    return dicom_dir, dont_unpack, scan_type

def run_info(info_file, dicom_dir, dont_unpack=[]):
    """Read an mri_parse_scdmdir infofile and determine run types."""
    from os.path import join as pjoin
    import numpy as np
    from unpack_dicoms import is_moco, get_fmap_type


    # Initialize the output lists
    files = []
    groups =[]
    seqs = []
    names = []
    ftypes = []
    args = []

    # Read the file into a numpy array
    seqfile = np.genfromtxt(info_file, dtype=object)

    # Go through the file and build the output lists
    for line in seqfile:
        seqn = int(line[2])

        # dont_unpack should be list of ints
        if seqn in dont_unpack:
            continue

        # Basic info
        dcmfile = line[1]
        name = line[12]
        slices, tps = tuple([int(dim) for dim in line[8:10]])
        
        # Set the extra arg to a null string
        # As of now it only gets used for FLASH
        arg = ""
        
        # Fill in output lists based on the acquisition
        if (slices==176) and (tps==1) and name.startswith("T1_MPRAGE"):
            # We want to convert the mprage twice (once to nii, once to mgh)
            # This is hacky, but let's update the output lists directly here
            # once, and then also set the variables get updated at the end of
            # this big if elif block
            files.append(dcmfile)
            groups.append("mri")
            seqs.append(seqn)
            names.append("001")
            ftypes.append("mgz")
            args.append("")
            
            # And this stuff gets added later 
            group = "structural"
            ftype = "niigz"
            name = "mprage"

        elif name.startswith("field_mapping"): 
            fmaptype = get_fmap_type(pjoin(dicom_dir, dcmfile))
            if name == "field_mapping":
                name = "func_%s_fm"%fmaptype
            else:
                name = "".join([name.replace("field_mapping_","").lower(), "_%s_fm"%fmaptype])
            group = "fieldmaps"
            ftype = "niigz"

        elif (tps==70) and name.startswith("DIFFUSION"):
            group = "dwi"
            name = "diffusion"
            ftype = "niigz"

        elif (slices==176) and name.startswith("gre_mgh_multiecho"):
            # Get the flip angle from the acquisition name
            alpha = int(name[18:-12])

            if tps == 1:
                # RMS image gets handled normally
                name = "flash_%02d-rms"%alpha
                group = "flash"
                ftype = "mgz"

            elif tps == 8:
                # For the echos, modify the output lists directly and then bail out
                for echo in range(8):
                    files.append(dcmfile)
                    seqs.append(seqn)
                    names.append("flash_%02d-%d"%(alpha, echo))
                    groups.append("flash")
                    ftypes.append("mgz")
                    args.append("-nth %d"%echo)
                continue 
            else:
                # Getting here would be really weird
                # Possibly exception worthy?
                continue

        elif (name=="ge_func_2x2x2_Resting") and (tps==62) and not is_moco(pjoin(dicom_dir, dcmfile)):
            group = "bold"
            name = "Resting"
            ftype = "niigz"

        elif ("ge_func" in name) and not is_moco(pjoin(dicom_dir, dcmfile)):
            # This seems like a hack?
            fullacqs = dict(NBack=198, IQ=244, MOT_Block=248, MOT_Jitter=302)
            full = False
            for acqname, trs in fullacqs.items():
                if name.startswith(acqname) and tps==trs:
                    full = True
            if not full:
                continue

            group = "bold"
            name = name
            ftype = "niigz"
            
        else:
            continue

        # If we make it here, we want this sequence
        # So update the output lists
        files.append(dcmfile)
        groups.append(group)
        seqs.append(seqn)
        names.append(name)
        ftypes.append(ftype)
        args.append(arg)

        # Now get rid of the sequence variables
        # I think this should help prevent silly bugs
        del dcmfile, group, seqn, name, ftype, arg
        
    return files, groups, seqs, names, ftypes, args

def is_moco(dcmfile):
    """Determine if a run has on-line motion correction."""
    import sys
    import subprocess
    cmd = ['mri_probedicom', '--i', dcmfile, '--t', '8', '103e']
    proc  = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    stdout = proc.communicate()[0]
    return stdout.strip().startswith('MoCoSeries')

def get_fmap_type(dcmfile):
    """Determine the type of a fieldmapping dicom."""
    import sys
    import subprocess
    cmd = ["mri_probedicom", "--i", dcmfile, "--t", "8", "8"]
    proc = subprocess.Popen(cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout = proc.communicate()[0]
    if stdout.strip()[-4] == "M":
        return "mag"
    elif stdout.strip()[-4] == "P":
        return "phase"
    else:
        raise TypeError("Could not determine fieldmap type of %s"%dcmfile)

def get_dicom_path(dicom_dir, dicom_files, args):
    """Given inputs, return a list of paths to dicom files."""
    from os.path import join
    out_files = []
    for i, f in enumerate(dicom_files):
        if args[i]:
            echo = int(args[i][-1])
            f = f.replace("-1.dcm","-%d.dcm"%(1+(echo*176)))
        out_files.append(join(dicom_dir, f))
    return out_files

def write_file_lists(dicom_dir, dicom_files, seq_numbers, groups, args):
    """Write text files containing a list of dicom files in each acquistion we want to convert."""
    from glob import glob
    from os import getcwd
    from os.path import join as pjoin

    fl_files = []
    
    for i, fname in enumerate(dicom_files):
        fhead = pjoin(dicom_dir, fname[:-5])
        if groups[i] == "flash" and args[i]:
            # Do multiecho flash files
            # This works because args[i] will be null string unless it's a FLASH echo
            echo = int(args[i][-1])
            flf = pjoin(getcwd(), "seq_%02d_echo_%d_files.txt"%(seq_numbers[i], echo))
            with open(flf, "w") as fobj:
                # Build the list manually
                start = 1 + (echo * 176)
                end = start + 176
                files = ["".join([fhead, str(f), ".dcm"]) for f in xrange(start, end)]
                fobj.write(" ".join(files))
            fl_files.append(flf)
        else:
            # Everything else gets written here
            flf = pjoin(getcwd(), "seq_%02d_files.txt"%seq_numbers[i])
            with open(flf, "w") as fobj:
                # Buld the list by a glob (maybe dangerous?)
                files = glob("".join([fhead, "*.dcm"]))
                fobj.write(" ".join(files))
            fl_files.append(flf)

    return fl_files

def rename_images(scan_type, converted_images, seq_numbers, ftypes, args):
    """Rename the list of images to reflect provenance info before datasinking."""
    from os import symlink, getcwd
    from os.path import join as pjoin

    extdict = dict(niigz="nii.gz", mgz="mgz")
    out_files = []
    for i, img in enumerate(converted_images):
        imgstring = "%s_seq_%02d"%(scan_type, seq_numbers[i])
        if args[i]:
            imgstring += "_echo_%d"%int(args[i][-1])
        dst = pjoin(getcwd(), ".".join([imgstring, extdict[ftypes[i]]]))
        symlink(img, dst)
        out_files.append(dst)
    return out_files

def rename_info_file(scan_type, info_file):
    from os import symlink, getcwd
    from os.path import join as pjoin

    dst = pjoin(getcwd(), "%s-dicominfo.txt"%scan_type)
    symlink(info_file, dst)
    return dst

def cache_run_info(data_dir, subject_list, types, dontunpack=[]):
    """Write a numpy array of run info to be used after dicoms are zipped."""
    for i, subj in enumerate(subject_list):
        type = types[i]
        dicom_dir = pjoin(data_dir, subj, "dicom", type)
        info_file = pjoin(data_dir, subj, "unpack", "%s-dicominfo.txt"%type)
        infoarr = np.array(run_info(info_file, dicom_dir, dontunpack))
        cache_file = pjoin(data_dir, subj, "unpack", "%s-infocache.npy"%type)
        np.save(cache_file, infoarr)

def load_info_cache(cache_file):
    """Load cached run info."""
    info = np.load(cache_file)
    return tuple([a.squeeze() for a in np.split(info,6)])

def make_symlinks(data_dir, subject_list, scan_types, dont_unpack=[]):
    """Create symlinks with heurtistic names to unpacked files."""

    for i, subj in enumerate(subject_list):
        scan_type = scan_types[i]
        subj_dir = pjoin(data_dir, subj)
        
        # Read the cached info data to get arrays of information
        cache_file = pjoin(data_dir, subj, "unpack", "%s-infocache.npy"%scan_type)
        files, groups, seqs, names, ftypes, args = load_info_cache(cache_file)

        # Make a functional hash to deal with multiple runs
        bold_hash = dict(IQ=1,NBack=1,MOT_Block=1,MOT_Jitter=1)

        # Dict for extensions
        type_dict = dict(niigz="nii.gz", mgz="mgz")

        # Iterate through the acquistions that got converted
        for j, seq in enumerate(seqs):
            group = groups[j]
            name = names[j]
            extra_arg = args[j]

            seq = int(seq)
            # Define the source image
            if group == "flash" and extra_arg:
                echo = int(extra_arg[-1])
                src = pjoin(subj_dir, "nifti", "%s_seq_%02d_echo_%d.mgz"%(scan_type, seq, echo))
            else:
                src = pjoin(subj_dir, "nifti", "%s_seq_%02d.%s"%(scan_type, seq, type_dict[ftypes[j]]))

            # Build the dst name based on info
            dst_dir = pjoin(subj_dir, group)
            if group == "mri":
                dst_dir = pjoin(dst_dir, "orig")
                dst =  "001.mgz"
            
            elif group == "structural":
                dst = "mprage.nii.gz"

            elif group == "dwi":
                dst = "diffusion.nii.gz"
            
            elif group == "fieldmaps":
                dst = ".".join([name, "nii.gz"])

            elif group == "flash":
                # Echo information is in flash name field
                dst = ".".join([name, "mgz"])

            elif group == "bold":
                
                if name == "Resting":
                    dst = "Resting.nii.gz"
                else:
                    if name.startswith("MOT"):
                        par = "_".join(name.split("_")[:2])
                    else:
                        par = "_".join(name.split("_")[:1])

                    dst = "%s_run%d.nii.gz"%(par, bold_hash[par])

                    bold_hash[par] += 1

            if os.path.exists(src):
                if not os.path.exists(dst_dir):
                    os.makedirs(dst_dir)
                
                dst = pjoin(dst_dir, dst)

                if not os.path.lexists(dst):
                    os.symlink(src, pjoin(dst_dir, dst))
                else:
                    print "Link %s already exists"%dst
            else:
                print "Source file %s not found"%src

            del group, name, extra_arg

def preprocess_anatomicals(data_dir, subject_list, types, stages):
    """Write and submit a Sun Grid Engine script to preprocess anatomical images."""

    for i, subj in enumerate(subject_list):
        scan_type = types[i]
        # Just bail out of this subject had a functional scan
        if scan_type == "func":
            continue
        proc_list = []
        subj_dir = pjoin(data_dir, subj)
        if "recon" in stages:
            # Check that the source image exists
            if not os.path.exists(pjoin(subj_dir, "mri", "orig", "001.mgz")):
                print "recon-all requested for %s, but MPRAGE image does not exist"%subj
                continue
            # Make sure a recon hasn't already been started for this subject
            # (We'll need to add functionality to allow for forcing a restart)
            if os.path.exists(pjoin(subj_dir, "scripts", "recon-all-status.log")):
                print "recon-all requested for %s, but recon-all.status.log already exists"%subj
                continue
            
            # Add the recon-all command to the processing list
            proc_list.append("recon-all -s %s -all"%subj)

        if "dwi" in stages:
            # Figure out where the diffusion DICOM is
            dicom_dir = pjoin(subj_dir, "dicom", scan_type)
            cache_file = pjoin(data_dir, subj, "unpack", "%s-infocache.npy"%scan_type)
            info = load_info_cache(cache_file)
            # Info pops out as a tuple
            names = info[3]
            dicom_files = info[0]
            try:
                dwi_index = [s for s, n in enumerate(names) if n == "diffusion"][0]
            except IndexError:
                print "Could not find diffusion DICOM in scan info file for %s"%subj
                continue
            
            # Here's the first element of the diffusion DICOM series
            dwi_dicom = pjoin(dicom_dir, dicom_files[dwi_index])

            # Specify the output dir
            dwi_dir = pjoin(subj_dir, "dwi")

            # Add the dt_recon command to the processing list
            proc_list.append("dt_recon --i %s --s %s --o %s --no-tal"%(dwi_dicom, subj, dwi_dir))

        if "norm" in stages:
            # Check to make sure the recon-all source image exists at least.
            # The actual source images for normalization actually get created
            # during the reconstruction process, but this should be a bit safer.
            if not os.path.exists(pjoin(subj_dir, "mri", "orig", "001.mgz")):
                print "Spatial normalization requested for %s, but MPRAGE image does not exist"%subj
                continue
           
            # Get a ref to the normalization script
            norm_script = "/mindhive/gablab/fluid/Nipype_Code/fluid_normalize.py"
            # Add the normalization command to the processing list
            proc_list.append("python %s -s %s"%(norm_script, subj))

        # Now, if we're supposed to be doing anything, tell the SGE to do it
        if proc_list:
            submit_to_sge(subj, proc_list)
            
            
def submit_to_sge(subject_id, proc_list):
    """Submit a list of commands to the Sun Grid Engine."""
    # Make sure a directory exists for SGE stuff
    sge_dir = pjoin(os.environ["HOME"], "sge")
    if not os.path.exists(sge_dir):
        os.mkdir(sge_dir)

    # Open a qsub script and fill it with processing info
    job_name = "%s_recon"%subject_id
    qsub_script = pjoin(sge_dir, "%s_%s.sh"%(job_name, time.time()))
    with open(qsub_script, "w") as q:
        # qsub arguments
        header = ["#! /bin/bash",
                  "#$ -V",
                  "#$ -cwd",
                  "#$ -N %s"%job_name,
                  "#$ -r y",
                  "#$ -S /bin/bash",
                  "\n"]
        q.write("\n".join(header)) 
        
        # Brain processing
        q.write("\n".join(proc_list))
    
    print "Submitting %s to Sun Grid Engine"%job_name
    qsub_cmd = " ".join(["cd", sge_dir, ";", "qsub", "-q", "long.q", qsub_script])

    # And actually submit the job
    proc = subprocess.Popen(qsub_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                            env=os.environ, shell=True, cwd=sge_dir)
    
    proc.communicate()

def compress_dicoms(data_dir, subject_list, types):
    """Tar and gzip DICOM directories for a list of subjects."""

    for i, subj in enumerate(subject_list):
        scan_type = types[i]
        
        # Move to the subject's directory 
        orig_dir = os.getcwd()
        dicom_dir = pjoin(data_dir, subj, "dicom")
        os.chdir(dicom_dir)

        print "Compressing dicom directory for %s"%subj

        # Tar the dicoms
        tar_file = "%s.tar"%scan_type
        dicom_glob = "/".join([scan_type, "*.dcm"])
        tar_cmd = "tar -cWf %s %s"%(tar_file, dicom_glob)
        os.system(tar_cmd)
        
        # Compress them
        zip_cmd = "gzip -f %s"%tar_file
        os.system(zip_cmd)

        # And, if things seemed to work, remove the directory of files
        zip_file =  "%s.gz"%tar_file
        if os.path.exists(zip_file):
            shutil.rmtree(scan_type)
        else:
            print "%s not created; DICOM directory will not be removed"%zip_file
        
        # Go back home
        os.chdir(orig_dir)

def parse_cmdline(arglist):
    """Parse a list of arguments."""
    parser = argparse.ArgumentParser(usage=("unpack_dicoms.py [options]"))
    parser.add_argument("-subjects", nargs="*", metavar="subjid",
                        help="list of subject ids to unpack")
    parser.add_argument("-procs", nargs="*", metavar="stages",
                        choices=["fetch","unpack","cache","link","anat","recon","dwi","norm","zip"],
                        help="processing stages to run")
    parser.add_argument("-type", metavar="scantype", choices=["func", "struct", "full"],
                        help="scan type - gets from # of dicoms if not specified")
    parser.add_argument("-dontunpack", nargs="*", metavar="num", default=[], type=int,
                        help="sequence number of run(s) to skip")
    parser.add_argument("-relink", action="store_true",
                        help="overwrite any old heuristic links")
    parser.add_argument("-workingdir", help="specify working directoy")
    parser.add_argument("-ipython", action="store_true",
                        help="run in parallel using IPython")
    
    if len(arglist) < 2:
        arglist.insert(0,"-h")
    args = parser.parse_args(arglist)

    return args

if __name__ == "__main__":
    main(sys.argv[1:])
