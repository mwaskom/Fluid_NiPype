#! /usr/bin/env python
"""
"""
import os
import sys
from os.path import join as pjoin
import argparse

from nipype.interfaces import io
import nipype.interfaces.utility as util   
import nipype.pipeline.engine as pe         

from workflows.restingstate import create_resting_workflow

import fluid_utility as flutil

def main(arglist):

    # Define paths
    project_dir = "/mindhive/gablab/fluid/"
    data_dir = pjoin(project_dir, "Data")
    analysis_dir = pjoin(project_dir, "Analysis/Nipype/resting")
    working_dir = pjoin(project_dir, "Analysis/Nipype/workingdir/")

    # Set the Freesurfer Subjects Directory 
    os.environ["SUBJECTS_DIR"] = data_dir

    # Smoothing kernel
    # Currently used for both volume and surface
    smooth_fwhm = 6 

    # Parse the argument list
    args = parse_args(arglist)

    # Figure out the subjects
    if args.subjects is not None:
        if os.path.isfile(args.subjects[0]):
            s_list = open(args.subjects[0]).read().strip().split()
        else:
            s_list = args.subjects
    else:
        s_list = []

    # Get the top-level workflow and in/out nodes
    wf, wf_input, wf_output = create_resting_workflow()
    
    # Set the top level working dir
    wf.base_dir = working_dir

    # Subject source node
    subjectsource = pe.Node(util.IdentityInterface(fields=["subject_id"]),
                            iterables=("subject_id", s_list),
                            name="subjectsource")
    
    # Say where to find the input data
    datasource = pe.Node(io.DataGrabber(infields=["subject_id"],
                                        outfields=["timeseries", "warpfield"],
                                        base_directory=data_dir),
                         name="datasource")
    datasource.inputs.template = "%s/%s/%s.nii.gz"
    datasource.inputs.template_args = dict(timeseries=[["subject_id", "bold", "Resting"]],
                                           warpfield=[["subject_id", "normalization", "warpfield"]])
    wf.connect(subjectsource, "subject_id", datasource, "subject_id")

    # Connect the data inputs
    flutil.connect_inputs(wf, datasource, wf_input)

    # Set the smoothing and subject_id inputs
    wf_input.inputs.smooth_fwhm = smooth_fwhm
    wf.connect(subjectsource, "subject_id", wf_input, "subject_id")

    # Set up the datasink
    datasink = pe.Node(io.DataSink(base_directory=analysis_dir,
                                   regexp_substitutions=[(r"_\w+\d", ""),
                                                         (r"_hemi_[lr]h", "")]),
                       name="datasink")
    flutil.subject_container(wf, subjectsource, datasink)
    flutil.sink_outputs(wf, wf_output, datasink, "@") 

    # Archive crashdumps
    flutil.archive_crashdumps(wf)

    # Figure out which engine to use
    plugin_args = dict()
    if args.ipython:
        plugin = "IPython"
    elif args.torque:
        plugin = "PBS"
        plugin_args['qsub_args'] = "-q gablab"
    elif args.linear:
        plugin = "Linear"
    else:
        plugin = "MultiProc"
        plugin_args['n_procs'] = args.nprocs

    # Process the data!
    wf.run(plugin=plugin, plugin_args=plugin_args)
        
def parse_args(arglist):

    parser = argparse.ArgumentParser(description="Resting-state processing for GFluid project")
    parser.add_argument("-subjects", nargs="*",
                        metavar="subject_id",
                        help="process subject(s)")
    parser.add_argument("-ipython",action="store_true",
                        help="run in parallel using IPython")
    parser.add_argument("-torque", action="store_true",
                    help="run in parallel using Torque")
    parser.add_argument("-linear",action="store_true",
                        help="run in series")
    parser.add_argument("-nprocs", type=int, default=4, 
                        help="number of MultiProcessing procs")
    args = parser.parse_args(arglist)

    return args

if __name__ == "__main__":
    main(sys.argv[1:])


