#! /usr/bin/env python
"""
Prepare fieldmap images.

Requires Gf unpacking pipeline to have been run.
Dwell time and phase-map echo deltas are hardcoded into this script for now.
Inspired by Doug Greve epidewarp.fsl script.
And by inspired by, I mean "is a Python port of."
"""
import os
import sys
import shutil
import argparse
import subprocess
from tempfile import mkdtemp

def main():

    # Hardcode the directory structure
    data_dir = "/mindhive/gablab/fluid/Data"

    # Fieldmapes are prefixed with longer names than we'll use here
    prefixdict = dict(rest="resting",func="func",dti="diffusion")

    # Hardcode the EPI dwell-time and fieldmap TE difference
    dwelldict = dict(rest=0.7,func=0.5,dti=0.8) 
    for subj in args.subjects:
        mapdir = os.path.join(data_dir, subj, "fieldmaps")
        if not os.path.exists(mapdir):
            print "Subject %s has no fieldmap directory"%subj
            continue
        for protocol in args.protocols:
            phase_image = prefixdict[protocol]+"_phase_fm"
            mag_image = prefixdict[protocol]+"_mag_fm"
            for img in phase_image, mag_image:
                if not os.path.exists(os.path.join(mapdir, img + ".nii.gz")):
                    print "Fieldmap image %s does not exist"%img
                    continue
            if args.verbose:
                print "Running prepmap function for %s protocol"%protocol
            prepmap(mapdir,phase_image,mag_image,protocol,dwelldict[protocol],2.46)

def prepmap(mapdir, phase_image, mag_image, prefix, dwell_time, te_diff):

    try:
        # Get situated
        origdir = os.getcwd()
        os.chdir(mapdir)

        # Create a temporary directory for stuff
        if args.tempdir is None:
            tmp = mkdtemp()
        else:
            tmp = args.tempdir
            if not os.path.exists(tmp):
                os.mkdir(tmp)

        # Skullstrip the mag image
        brain = os.path.join(tmp, "brain")
        brainmask = os.path.join(tmp, "brain_mask")
        runcmd(["bet", mag_image, brain, "-m"])

        # Dilate the resulting brainmask to get a headmask
        headmask = os.path.join(tmp, "head_mask")
        runcmd(["fslmaths", brainmask, "-dilM -dilM -dilM", headmask])

        # Rescale the phase-difference image to -pi<vox<pi
        scaledphase = os.path.join(tmp, "scaled_phase")
        runcmd(["fslmaths -dt float", phase_image, "-sub 2047.5 -mul 0.00153436", 
                scaledphase, "-odt float"])

        # Unwrap the phase image using prelude
        unwrappedphase = os.path.join(tmp, "unwrapped_phase")
        runcmd(["prelude -p", scaledphase, "-a", mag_image, "-o", unwrappedphase, "-f -m", headmask])

        # Create a dummy phase echo and merge it with the phase diff map
        dummy = os.path.join(tmp, "dummy_phase")
        runcmd(["fslmaths", unwrappedphase, "-mul 0", dummy])
        mergedphase = os.path.join(tmp, "merged_phase")
        runcmd(["fslmerge -t", mergedphase, dummy, unwrappedphase])

        # Create the voxel-shift map
        origvsm = os.path.join(tmp, "orig_vsm")
        magdw = os.path.join(tmp, "mag_dw")
        runcmd(["fugue -i", mag_image, "-u", magdw, "-p", mergedphase, 
                "--dwell=%.6f"%(dwell_time*.001), "--asym=%.6f"%(te_diff*.001),
                "--mask=%s"%brainmask, "--saveshift=%s"%origvsm,
                "--smooth2=2"])

        # Demean the voxel-shift map
        vsmmean = runcmd(["fslstats", origvsm, "-k", brainmask, "-m"], return_stdout=True)
        vsm = prefix + "_voxel_shift_map"
        runcmd(["fslmaths", origvsm, "-sub", vsmmean, "-mul", brainmask, vsm])

        # Forard warp the skullstripped magnitude volume
        warpedmag = prefix+"_warped_mag"
        runcmd(["fugue -i", brain, "-w", warpedmag, "--loadshift=%s"%vsm, "--mask=%s"%brainmask])
    
    finally:
        # Clean-up
        os.chdir(origdir)
        if args.cleanup:
            shutil.rmtree(tmp)
        else:
            print "Not removing tempdir at %s"%tmp

def runcmd(cmdline, return_stdout=False):
    """Use the subprocess.Popen class to run a command"""
    if isinstance(cmdline, list): 
        cmd = " ".join(cmdline)
    else:
        cmd = cmdline
    if args.verbose:
        print cmd

    proc = subprocess.Popen(cmd, 
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            env=os.environ,
                            shell=True,
                            cwd=os.getcwd())
               

    stdout, stderr = proc.communicate()

    if stderr:
        raise RuntimeError(stderr)

    if args.verbose:
        print stdout
    if return_stdout:
        return stdout.strip()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-subjects", nargs="*", help="subjects (requires Gf unpacking)")
    parser.add_argument("-protocols", nargs="*", help="func, dti, or rest")
    parser.add_argument("-tempdir", help="use specified directory and don't clean up")
    parser.add_argument("-nocleanup", action="store_false", dest="cleanup",
                        help="do not remove temp dir")
    parser.add_argument("-verbose", action="store_true", help="verbose output")
    if len(sys.argv) == 1:
        sys.argv.insert(1,"-h")
    args = parser.parse_args()
    if args.protocols==["all"]:
        args.protocols = ["func","dti","rest"]
    if args.tempdir is not None:
        args.cleanup = False
    main()
