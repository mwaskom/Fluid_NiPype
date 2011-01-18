#! /usr/bin/env python
"""
Main interface for GFluid Freesurfer Morphometry Nipype code.
"""
import os
from os.path import join as pjoin
import sys
import argparse

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util

# Parse command line arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-subjects", nargs="*",
                    metavar="subject_id",
                    help="process subject(s)")
parser.add_argument("-workflows", nargs="*",
                    metavar="WF",
                    help="which workflows to run")
parser.add_argument("-fwhm", nargs="*",
                    type=float,
                    help="smoothing kernel")
parser.add_argument("-metric", nargs="*",
                    help="metric to use")
parser.add_argument("-inseries",action="store_true",
                    help="force running in series")
args = parser.parse_args()

for metric in args.metric:
    if metric not in ["thickness", "jacobian_white", "curv", "sulc"]:
        sys.exit("%s is not a valid metric"%metric)


project_dir = "/mindhive/gablab/fluid/"
data_dir = pjoin(project_dir, "Data")
analysis_dir = pjoin(project_dir, "Analysis/Nipype/morphometry/")
working_dir = pjoin(project_dir, "Analysis/Nipype/workingdir", "morphometry")

try:
    os.mkdir(analysis_dir)
except OSError:
    pass

subjsource = pe.Node(util.IdentityInterface(fields=["subject_id"]),
                     iterables=("subject_id", args.subjects),
                     name="subjectsource")

hemisource = pe.Node(util.IdentityInterface(fields=["hemi"]),
                     iterables=("hemi", ["lh", "rh"]),
                     name="hemisource")

fwhmsource = pe.Node(util.IdentityInterface(fields=["fwhm"]),
                     iterables=("fwhm", args.fwhm),
                     name="fwhmsource")

metricsource = pe.Node(util.IdentityInterface(fields=["metric"]),
                       iterables=("metric",args.metric),
                       name="metricsource")

surfgrabber = pe.Node(nio.DataGrabber(infields=["subject_id", "metric", "hemi"],
                                      outfield=["metric_map"],
                                      base_directory=data_dir,
                                      template="%s/surf/%s.%s"),
                      name="surfgrabber")
surfgrabber.inputs.template_args = dict(metric_map=[["subject_id", "hemi", "metric"]])

transform = pe.Node(fs.SurfaceTransform(target_subject="fsaverage",
                                        target_type="mgz"),
                        name="surftransform")

smooth = pe.Node(fs.SurfaceSmooth(subject_id="fsaverage"),
                 name="smooth")

surfsink = pe.Node(nio.DataSink(base_directory=analysis_dir,
                                parameterization=False),
                   name="surfsink")

preproc = pe.Workflow(name="preproc", base_dir=working_dir)

preproc.connect([
    (subjsource,    surfgrabber,  [("subject_id", "subject_id")]),
    (hemisource,    surfgrabber,  [("hemi", "hemi")]),
    (metricsource,  surfgrabber,  [("metric", "metric")]),
    (surfgrabber,   transform,    [("metric_map", "source_file")]),
    (subjsource,    transform,    [("subject_id", "source_subject")]),
    (hemisource,    transform,    [("hemi", "hemi")]),
    (transform,     smooth,       [("out_file", "in_file")]),
    (fwhmsource,    smooth,       [("fwhm", "fwhm")]),
    (hemisource,    smooth,       [("hemi", "hemi")]),
    (smooth,        surfsink,     [("out_file", "@surf")]),
    (subjsource,    surfsink,     [("subject_id", "container")]),
    ])

mergegrabber = pe.Node(nio.DataGrabber(infields=["hemi", "metric", "fwhm"],
                                       outfield=["metric_maps"],
                                       sort_filelist=True,
                                       base_directory=analysis_dir,
                                       template="gf??/%s.%s.fsaverage_smooth%d.mgz"),
                      name="mergegrabber")
mergegrabber.inputs.template_args = dict(metric_maps=[["hemi", "metric", "fwhm"]])
mergegrabber.overwrite=True

merger = pe.Node(fs.Concatenate(),name="merger")

rename = pe.Node(util.Merge(3), name="rename")

mergesink = pe.Node(nio.DataSink(base_directory=analysis_dir,
                                 parameterization=False),
                    name="mergesink")

merge = pe.Workflow(name="merge", base_dir=working_dir)

def get_formatted_name(namelist):
    
    return [("concat_output", "%s.%s.fsaverage.fwhm%d"%tuple(namelist))]

merge.connect([
    (hemisource,    mergegrabber,  [("hemi", "hemi")]),
    (metricsource,  mergegrabber,  [("metric", "metric")]),
    (fwhmsource,    mergegrabber,  [("fwhm", "fwhm")]),
    (mergegrabber,  merger,        [("metric_maps", "in_files")]),
    (merger,        mergesink,     [("concatenated_file", "group.@surface")]),
    (hemisource,    rename,        [("hemi", "in1")]),
    (metricsource,  rename,        [("metric", "in2")]),
    (fwhmsource,    rename,        [("fwhm", "in3")]),
    (rename,        mergesink,     [(("out", get_formatted_name), "substitutions")])
    ])

if __name__ == "__main__":
    if "preproc" in args.workflows:
        preproc.run(inseries=args.inseries)
    if "merge" in args.workflows:
        merge.run(inseries=args.inseries)
