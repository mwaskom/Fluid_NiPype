#! /usr/bin/env python
import os
import argparse
from glob import glob
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe

from workflows.vbm_preproc import preproc

import fluid_utility as flutil

# Parse command line arguments
parser = argparse.ArgumentParser(description="VBM-style DTI processing for GFluid project")
parser.add_argument("-subjects", nargs="*",
                    metavar="subject_id",
                    help="process subject(s)")
parser.add_argument("-workflows", nargs="*",
                    metavar="WF",
                    help="which workflows to run")
parser.add_argument("-norun",dest="run",action="store_false",
                    help="do not run the pypeline")
parser.add_argument("-inseries",action="store_true",
                    help="force running in series")
args = parser.parse_args()

project_dir = "/mindhive/gablab/fluid/"
data_dir = "/mindhive/gablab/fluid/Data"
analysis_dir = "/mindhive/gablab/fluid/Analysis/Nipype/dti_vbm"
working_dir = "/mindhive/gablab/fluid/Analysis/Nipype/workingdir/dti_vbm"
report_dir = "/mindhive/gablab/fluid/Analysis/Nipype/dti_vbm/report"

# Get all of the pre subjects as a default subject list 
if args.subjects is None:
    paths = glob(os.path.join(data_dir, "gf??"))
    paths = [p for p in paths if os.path.exists(os.path.join(p, "bold", "Resting.nii.gz"))]
    subject_list = [p.split(os.path.sep)[-1] for p in paths]
else:
    subject_list = args.subjects

print "Subjects: ", " ".join(subject_list)

if args.workflows is None:
    args.workflows = []
elif args.workflows == ["all"]:
    args.workflows = ["preproc"]

# Subject source node
# -------------------
subjectsource = pe.Node(util.IdentityInterface(fields=["subject_id"]),
                        iterables=("subject_id", subject_list),
                        name="subjectsource")

# Preprocessing
# -------------

# Get preproc input and output
preproc_input = preproc.get_node("inputspec")
preproc_output = preproc.get_node("outputspec")
preproc_report = preproc.get_node("report")

# Define a data source to grab input data
preproc_source = pe.Node(nio.DataGrabber(infields=["subject_id"],
                                        outfields=["stat_image", "affine_mat", "warpfield"],
                                        base_directory=data_dir),
                         name="preprocsource")

preproc_source.inputs.template = "%s/%s/%s"
preproc_source.inputs.template_args = dict(stat_image = [["subject_id", "dwi", "fa.nii"]],
                                           affine_mat = [["subject_id", "dwi", "fsl_reg.mat"]],
                                           warpfield = [["subject_id", "normalization", "warpfield.nii.gz"]])

                    
# Set the smoothing value
preproc_input.inputs.smooth_fwhm = 10

# Connect the inputs
preproc.connect(subjectsource, "subject_id", preproc_source, "subject_id")
flutil.connect_inputs(preproc, preproc_source, preproc_input)

# Sink the outputs
preprocsink = pe.Node(nio.DataSink(base_directory=analysis_dir), name="datasink")
preprocreport = pe.Node(nio.DataSink(base_directory=report_dir), name="reportsink")

preprocsinksubnode = pe.Node(util.Merge(len(preproc_output.outputs.__dict__.keys())),
                             name="sinksubnode")
preprocreportsubnode = pe.Node(util.Merge(len(preproc_report.outputs.__dict__.keys())), 
                               name="reportsubnode")

# Connect the outputs to the sink
flutil.subject_container(preproc, subjectsource, preprocsink)
flutil.subject_container(preproc, subjectsource, preprocreport)
flutil.sink_outputs(preproc, preproc_output, preprocsink, "preproc")
flutil.sink_outputs(preproc, preproc_report, preprocreport, "preproc")


flutil.get_output_substitutions(preproc, preproc_output, preprocsinksubnode)
flutil.get_output_substitutions(preproc, preproc_report, preprocreportsubnode)

# Make substitutions 

# Change the final image name to something more informative
preprocsinknodesubs = [("stat_image.nii.gz", "fa.nii.gz")]

flutil.set_substitutions(preproc, preprocsink, preprocsinksubnode, preprocsinknodesubs)
flutil.set_substitutions(preproc, preprocreport, preprocreportsubnode, [])

# Set the base directory
preproc.base_dir = os.path.join(working_dir)

# Archive crashdumps
flutil.archive_crashdumps(preproc)


# Run the pipeline
# ----------------

if __name__ == "__main__":
    if args.run:
        if __file__ == "fluid_resting.py":
            report_script = "python /mindhive/gablab/fluid/Nipype_Code/reporting/build_report.py "
        else:
            report_script = "python /dev/null "
        
        if "preproc" in args.workflows:
            preproc.run(inseries=args.inseries)
            os.system(report_script + " ".join(subject_list))
        
        if "reg" in args.workflows:
            registration.run(inseries=args.inseries)
            os.system(report_script + " ".join(subject_list))
