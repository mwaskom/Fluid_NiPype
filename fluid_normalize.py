#! /usr/bin/env python

import os
import argparse

from nipype.interfaces import io
from nipype.interfaces import fsl
from nipype.interfaces import freesurfer as fs
from nipype.interfaces import utility as util
from nipype.interfaces import base
from nipype.interfaces.fsl import base as fslbase
from nipype.utils.filemanip import fname_presuffix
from nipype.pipeline import engine as pe

from fluid_utility import archive_crashdumps

# Parse command line arguments
# ----------------------------
parser = argparse.ArgumentParser(description="Perform spatial normalization to the MNI template.")

parser.add_argument("-subjects", nargs="*",
                    metavar="subject_id",
                    help="process subject(s)")
parser.add_argument("-nopype", dest="pype", action="store_false",
                    help="don't run the normalization workflow")
parser.add_argument("-parallel", dest="inseries", action="store_false", 
                    help="run workflows in parallel if cluster engine is availible")
args = parser.parse_args()

# Define some paths
# -----------------
project_dir = "/mindhive/gablab/fluid"
data_dir = os.path.join(project_dir, "Data")
working_dir = os.path.join(project_dir, "Analysis/Nipype/workingdir/normalization")

# Get target images
target_brain = fsl.Info.standard_image("avg152T1_brain.nii.gz")
target_head =  fsl.Info.standard_image("avg152T1.nii.gz")
target_mask = fsl.Info.standard_image("MNI152_T1_2mm_brain_mask_dil.nii.gz")
fnirt_cfg = os.path.join(os.environ["FSLDIR"], "etc/flirtsch/T1_2_MNI152_2mm.cnf")
mni_orient = ("RL", "PA", "IS")

# Define a CheckReg interface here for now (will likely migrate to Nipype source)
# -------------------------------------------------------------------------------
class CheckRegInput(fslbase.FSLCommandInputSpec):
    
    in_file = base.File(exists=True, argstr="%s", position=1)
    out_file = base.File(genfile=True, argstr="%s", position=2)

class CheckRegOutput(base.TraitedSpec):

    out_file = base.File(exists=True)

class CheckReg(fslbase.FSLCommand):

    _cmd = "check_mni_reg"
    input_spec = CheckRegInput
    output_spec = CheckRegOutput

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["out_file"] = fname_presuffix(self.inputs.in_file,
                                              suffix="_to_mni.png",
                                              use_ext=False,
                                              newpath=os.getcwd())
        return outputs

    def _gen_filename(self, name):
        if name == "out_file":
            return self._list_outputs()[name]
        return None

# Set up the workflow
# -------------------

# Subject source node
subjectsource = pe.Node(util.IdentityInterface(fields=["subject_id"]),
                        iterables = ("subject_id", args.subjects),
                        name = "subjectsource")

# Grab recon-all outputs
datasource = pe.Node(io.DataGrabber(infields=["subject_id"],
                                    outfields=["brain", "head"],
                                    base_directory=data_dir,
                                    template="%s/mri/%s.mgz"),
                     name="datagrabber")
datasource.inputs.template_args = dict(brain=[["subject_id", "norm"]],
                                       head=[["subject_id", "nu"]])

# Convert images to nifti storage and float representation
cvtbrain = pe.Node(fs.MRIConvert(out_type="niigz", out_datatype="float"),
                   name="convertbrain")

cvthead = pe.Node(fs.MRIConvert(out_type="niigz", out_datatype="float"),
                  name="converthead")

# FLIRT brain to MNI152_brain
flirt = pe.Node(fsl.FLIRT(reference=target_brain),
                name="flirt")
sw = [-180, 180]
for dim in ["x","y","z"]:
    setattr(flirt.inputs, "searchr_%s"%dim, sw)

# FNIRT head to MNI152
fnirt = pe.Node(fsl.FNIRT(ref_file=target_head,
                          refmask_file=target_mask,
                          config_file=fnirt_cfg,
                          fieldcoeff_file=True),
                name="fnirt")

# Warp the images
warpbrain = pe.Node(fsl.ApplyWarp(ref_file=target_head,
                                  interp="spline"),
                    name="warpbrain")

warphead = pe.Node(fsl.ApplyWarp(ref_file=target_head,
                                 interp="spline"),
                    name="warphead")

# Generate a png summarizing the registration
checkreg = pe.Node(CheckReg(),
                   name="checkreg")

# Save relevant files to the data directory
datasink = pe.Node(io.DataSink(base_directory=data_dir,
                               parameterization = False,
                               substitutions=[("norm_", "brain_"),
                                              ("nu_", "T1_"),
                                              ("_out", ""),
                                              ("T1_fieldwarp", "warpfield"),
                                              ("brain_flirt.mat", "affine.mat")]),
                   name="datasink")

# Define and connect the workflow
# -------------------------------

normalize = pe.Workflow(name="normalize", 
                        base_dir=working_dir)

archive_crashdumps(normalize)

normalize.connect([
    (subjectsource,   datasource,   [("subject_id", "subject_id")]),
    (datasource,      cvtbrain,     [("brain", "in_file")]),
    (datasource,      cvthead,      [("head", "in_file")]),
    (cvtbrain,        flirt,        [("out_file", "in_file")]),
    (flirt,           fnirt,        [("out_matrix_file", "affine_file")]),
    (cvthead,         fnirt,        [("out_file", "in_file")]),
    (cvtbrain,        warpbrain,    [("out_file", "in_file")]),
    (fnirt,           warpbrain,    [("fieldcoeff_file", "field_file")]),
    (cvthead,         warphead,     [("out_file", "in_file")]),
    (fnirt,           warphead,     [("fieldcoeff_file", "field_file")]),
    (warpbrain,       checkreg,     [("out_file", "in_file")]),
    (subjectsource,   datasink,     [("subject_id", "container")]),
    (cvtbrain,        datasink,     [("out_file", "normalization.@brain")]),
    (cvthead,         datasink,     [("out_file", "normalization.@t1")]),
    (flirt,           datasink,     [("out_file", "normalization.@brain_flirted")]),
    (flirt,           datasink,     [("out_matrix_file", "normalization.@affine")]),
    (warphead,        datasink,     [("out_file", "normalization.@t1_warped")]),
    (warpbrain,       datasink,     [("out_file", "normalization.@brain_warped")]),
    (fnirt,           datasink,     [("fieldcoeff_file", "normalization.@warpfield")]),
    (checkreg,        datasink,     [("out_file", "normalization.@reg_png")]),
    ])

if __name__ == "__main__" and args.pype:
    normalize.run(inseries=args.inseries)
