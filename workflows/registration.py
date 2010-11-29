import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.fsl as fsl
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util


registration = pe.Workflow(name="registration")

# Define the inputs for the registation workflow
inputnode = pe.Node(util.IdentityInterface(fields=["subject_id",
                                                   "warpfield",
                                                   "mean_func",
                                                   "example_func",
                                                   "functional_mask",
                                                   "vol_timeseries",
                                                   "surf_timeseries"]),
                    name="inputspec")

func2anat = pe.MapNode(fs.BBRegister(contrast_type="t2",
                                     init="fsl",
                                     epi_mask=True,
                                     registered_file=True,
                                     out_fsl_file=True),
                       iterfield=["source_file"],
                       name="func2anat")

func2anatpng = pe.MapNode(fsl.Slicer(middle_slices=True,
                                     show_orientation=False,
                                     label_slices=False),
                          iterfield=["in_file"],
                          name="func2anatpng")

fssource  = pe.Node(nio.FreeSurferSource(subjects_dir=fs.Info.subjectsdir()),
                    name="fssource")

convert = pe.Node(fs.MRIConvert(out_type="niigz"), name="convertbrain")

mni152 = fsl.Info.standard_image("avg152T1.nii.gz")
mni152_brain = fsl.Info.standard_image("avg152T1_brain.nii.gz")

warpts = pe.MapNode(fsl.ApplyWarp(ref_file=mni152),
                    iterfield=["in_file", "premat"],
                    name="warptimeseries")

warpex = pe.MapNode(fsl.ApplyWarp(ref_file=mni152),
                    iterfield=["in_file", "premat"],
                    name="warpexample")

warpmask = pe.MapNode(fsl.ApplyWarp(ref_file=mni152),
                      iterfield=["in_file", "premat"],
                      name="warpmask")

warpmaskpng = pe.MapNode(fsl.Slicer(middle_slices=True,
                                    show_orientation=False,
                                    label_slices=False,
                                    image_edges=mni152_brain),
                         iterfield=["in_file"],
                         name="warpmaskpng")

meanwarp = pe.MapNode(fsl.ImageMaths(op_string="-Tmean", suffix="_mean"),
                      iterfield=["in_file"],
                      name="meanwarp")

exwarppng = pe.MapNode(fsl.Slicer(middle_slices=True,
                                  show_orientation=False,
                                  label_slices=False,
                                  image_edges=mni152_brain),
                       iterfield=["in_file"],
                       name="exwarppng")

outputnode = pe.Node(util.IdentityInterface(fields=["warped_timeseries",
                                                    "warped_example_func",
                                                    "warped_functional_mask",
                                                    "warped_mean_func",
                                                    "register"]),
                     name="outputspec")

report = pe.Node(util.IdentityInterface(fields=["func2anat", 
                                                "func2anat_cost", 
                                                "example_func_warp",
                                                "functional_mask_warp"])
                                                ]),
                 name="report")

registration.connect([
    (inputnode,    func2anat,      [("subject_id", "subject_id"),
                                    ("mean_func", "source_file")]),
    (inputnode,    warpts,         [("warpfield", "field_file"),
                                    ("vol_timeseries", "in_file")]),
    (inputnode,    warpex,         [("warpfield", "field_file"),
                                    ("example_func", "in_file")]),
    (inputnode,    warpmask,       [("warpfield", "field_file"),
                                    ("functional_mask", "in_file")]),
    (inputnode,    fssource,       [("subject_id", "subject_id")]),
    (fssource,     convert,        [("brain", "in_file")]),
    (convert,      func2anatpng,   [("out_file", "image_edges")]),
    (func2anat,    func2anatpng,   [("registered_file", "in_file")]),
    (func2anat,    warpts,         [("out_fsl_file", "premat")]),
    (func2anat,    warpex,         [("out_fsl_file", "premat")]),
    (func2anat,    warpmask,       [("out_fsl_file", "premat")]),
    (warpts,       meanwarp,       [("out_file", "in_file")]),
    (warpex,       exwarppng,      [("out_file", "in_file")]),
    (warpmask,     warpmaskpng,    [("out_file", "in_file")]),	
    (warpex,       outputnode,     [("out_file", "warped_example_func")]),
    (func2anat,    outputnode,     [("out_reg_file", "register")]),
    (warpts,       outputnode,     [("out_file", "warped_timeseries")]),
    (warpmask,     outputnode,     [("out_file", "warped_functional_mask")]),
    (meanwarp,     outputnode,     [("out_file", "warped_mean_func")]),
    (func2anat,    report,         [("min_cost_file", "func2anat_cost")]),
    (func2anatpng, report,         [("out_file", "func2anat")]),
    (warpmaskpng,  report,         [("out_file", "functionak_mask_warp"]),
    (exwarppng,    report,         [("out_file", "example_func_warp")])
    ])

