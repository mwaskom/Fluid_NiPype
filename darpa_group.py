import os
import re
import sys
import shutil
from copy import deepcopy
from datetime import datetime

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.spm as spm
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util

import argparse

if __name__ == "__main__":
    """ Handle command line arguments to control the analysis """
    parser = argparse.ArgumentParser(description="Main interface for GFluid NiPype code.")
    parser.add_argument("-paradigm", metavar="paradigm", required=True,
                        help="experimental paradigm")
    parser.add_argument("-subject", dest="subjects", 
                        metavar="subject_id", action="append",
                        help="run pypeline for subject(s)")
    parser.add_argument("-norun",dest="run",action="store_false",
                        help="do not run the pypeline")
    parser.add_argument("-nograph",dest="write_graph",action="store_false",
                        help="do not write a graph")
    parser.add_argument("-inseries",action="store_true",
                        help="force running in series")
    args = parser.parse_args()
    """ Dynamically import the experiment file """
else:
    class Foo(): pass
    args = Foo()
    args.paradigm = "nback"
    args.run = False
    args.write_graph = False
exp = __import__("%s_experiment" % args.paradigm)


gpipe = pe.Workflow(name="group_analysis")
working_output = os.path.join(
    os.path.abspath("../working"), args.paradigm)
gpipe.base_dir = working_output

copesource = pe.Node(util.IdentityInterface(fields=["contrast"]),
                     name="copesource",
                     iterables=("contrast", range(len(exp.contrasts))))

copegrabber = pe.Node(nio.DataGrabber(infields=["contrast"],outfields=["copes"]),
                      name="copegrabber")

copegrabber.inputs.base_directory = (
    "/mindhive/gablab/fluid/Analysis/DARPA_Visit/working/%s/level1/vol_fixed_fx/_fwhm_5.0/"%args.paradigm)
if args.paradigm == "nback":
    copegrabber.inputs.template = (
        "_subject_id_gf??/_index_0/flameo/mapflow/_flameo%d/stats/cope1.nii.gz")
else:
    copegrabber.inputs.template = (
        "_subject_id_gf??/flameo/mapflow/_flameo%d/stats/cope1.nii.gz")
    
copegrabber.inputs.template_args = dict(copes=[["contrast"]])

copemerge = pe.Node(fsl.Merge(dimension="t"), name="copemerge")

design = pe.Node(fsl.L2Model(), name="design")

mni152 = fsl.Info.standard_image("avg152T1_brain.nii.gz")

flamemask = pe.Node(fs.Binarize(min=10, in_file=mni152), name="flamemask")

flame = pe.Node(fsl.FLAMEO(run_mode="ols"), name="flame")

mni_brain = fsl.Info.standard_image("avg152T1.nii.gz")
overlayflame = pe.Node(interface=fsl.Overlay(stat_thresh=(2.5, 5),
                                                auto_thresh_bg=True,
                                                show_negative_stats=True,
                                                background_image=mni_brain),
                          name='overlayflame')

sliceflame = pe.Node(interface=fsl.Slicer(all_axial=True,
                                             image_width=750),
                        name='sliceflame')

hemisource = pe.Node(util.IdentityInterface(fields=["hemi"]),
                     iterables=("hemi",["lh","rh"]),
                     name="hemisource")

surfshots = pe.Node(fs.SurfaceScreenshots(subject="fluid_fsaverage",
                                          surface="semi-5",
                                          show_color_scale=True,
                                          show_gray_curv=True,
                                          six_images=True,
                                          overlay_range=(2.5,5),
                                          mni152_reg=True,
                                          tcl_script="/u2/mwaskom/darpa_screenshots.tcl"),
                     name="surfshots")

gpipe.connect([
    (copesource, copegrabber, [("contrast", "contrast")]),
    (copegrabber, copemerge,[("copes","in_files")]),
    (copegrabber, design, [(("copes", lambda x: len(x)), "num_copes")]),
    (copemerge,flame,[('merged_file','cope_file')]),
    (design,flame, [('design_mat','design_file'),
                    ('design_con','t_con_file'),
                    ('design_grp','cov_split_file')]),
    (flamemask,flame, [("binary_file", "mask_file")]),
    (flame,surfshots, [("zstats", "overlay")]),
    (hemisource, surfshots, [("hemi", "hemi")]),
    (flame,overlayflame,[('zstats','stat_image')]),
    (overlayflame,sliceflame,[('out_file', 'in_file')])
    ])

gpipe.config = dict(crashdump_dir="crashdumps/gpipe")

if args.run:
    gpipe.run(inseries=args.inseries)
if args.write_graph:
    gpipe.write_graph(graph2use="flat")
