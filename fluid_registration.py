import nipype.pipeline.engine as pe
import nippye.interfaces.io as nio
import nipype.interfaces.fsl as fsl
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util


registration = pe.Workflow(name="registration")

# Define the inputs for the registation workflow
inputnode = pe.Node(interface=util.IdentityInterface(fields=["subject_id",
                                                             "warpfield",
                                                             "example_func",
                                                             "vol_timeseries",
                                                             "surf_timeseries"]),
                    name="inputspec")

func2anat = pe.MapNode(interface=fs.BBRegister(contrast_type="t2",
                                               init="fsl",
                                               epi_mask=True,
                                               registered_file=True,
                                               out_fsl_file=True),
                       iterfield=["source_file"],
                       name="func2anat")

func2anatpng = pe.MapNode(interface=fsl.Slicer(middle_slices=True,
                                               show_orientation=False,
                                               label_slices=False),
                          iterfield=["in_file"],
                          name="func2anatpng")

fssource  = pe.Node(interface=nio.FreeSurferSource(subject_dir=fs.Info.subjectsdir()),
                    name="fssource")

mni152 = fsl.Info.standard_image("avg152T1.nii.gz")
mni152_brain = fsl.Info.standard_image("avg152T1_brain.nii.gz")

warpts = pe.MapNode(interface=fsl.ApplyWarp(ref_file=mni152),
                    iterfield=["in_file", "premat"],
                    name="warptimeseries")

warpex = pe.MapNode(interface=fsl.ApplyWarp(ref_file=mni152),
                    iterfield=["in_file", "premat"],
                    name="warpexample")

exfuncwarppng = pe.MapNode(interface=fsl.Slicer(middle_slices=True,
                                                show_orientation=False,
                                                label_slices=False,
                                                image_edges=mni152_brain),
                           iterfield=["in_file"],
                           name="exfuncwarppng")

outputnode = pe.Node(util.IdentityInterface(fields=["warped_timeseries",
                                                    "register"]),
                     name="outputspec")

report = pe.Node(util.IdentityInterface(fields=["func2anat", 
                                                "func2anat_cost", 
                                                "example_func_warp"
                                                ]),
                 name="report")

registration.connect([
    (inputnode,    func2anat,      [("subject_id", "subject_id"),
                                    ("example_func", "source_file")]),
    (inputnode,    warpts,         [("warpfield", "field_file"),
                                    ("vol_timeseries", "in_file")]),
    (inputnode,    warpex,         [("warpfield", "field_file"),
                                    ("example_func", "in_file")]),
    (fssource,     func2anatpng,   [("brain", "image_edges")]),
    (func2anat,    func2anatpng,   [("registered_file", "in_file")]),
    (func2anat,    warpts,         [("out_fsl_file", "premat")]),
    (func2anat,    warpex,         [("out_fsl_file", "premat")]),
    (warpex,       exfuncwarppng,  [("out_file", "in_file")]),
    (func2anat,    outputnode,     [("out_reg_file", "register")]),
    (warpts,       outputnode,     [("out_file", "warped_timeseries")]),
    (func2anat,    report,         [("min_cost_file", "func2anat_cost")]),
    (func2anatpng, report,         [("out_file", "func2anat")]),
    (exfuncwarppng,report,         [("out_file", "example_func_warp")])
    ])
