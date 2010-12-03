import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe


# Inputs
inputnode = pe.Node(util.IdentityInterface(fields=["cope", "varcope"]),
                    name="inputspec")

# Concatenate the Cope for each run
copemerge = pe.Node(fsl.Merge(dimension="t"),
                    name="copemerge")

# Concatenate the Varcope for each run
varcopemerge = pe.Node(fsl.Merge(dimension="t"),
                       name="varcopemerge")

# Set up a FLAMEO model
level2model = pe.Node(fsl.L2Model(),
                      name="l2model")

# Will probably want to change this so we can run in native space
brain_mask = fsl.Info.standard_image("MNI152_T1_2mm_brain_mask_dil.nii.gz")

# Run a fixed effects analysis in FLAMEO
flameo = pe.Node(fsl.FLAMEO(run_mode="fe", mask_file=brain_mask),
                 name="flameo")


# Again, will probably want to change this
mni_brain = fsl.Info.standard_image("avg152T1.nii.gz")

# Overlay the stats onto a background image
overlayflame = pe.Node(fsl.Overlay(stat_thresh=(2.3, 10),
                                                auto_thresh_bg=True,
                                                show_negative_stats=True,
                                                background_image=mni_brain),
                          name="overlayflame")

# Slice the overlaid statistical images
sliceflame = pe.Node(fsl.Slicer(image_width=872),
                        name="sliceflame")
sliceflame.inputs.sample_axial = 2

# Outputs
outputnode = pe.Node(util.IdentityInterface(fields=["stats"]), name="outputspec")

# Report
report = pe.Node(util.IdentityInterface(fields=["zstat"]), name="report")

fixed_fx = pe.Workflow(name="fixed_fx")

fixed_fx.connect([
    (inputnode,    copemerge,     [("cope", "in_files")]),
    (inputnode,    varcopemerge,  [("varcope", "in_files")]),
    (copemerge,    flameo,        [("merged_file","cope_file")]),
    (varcopemerge, flameo,        [("merged_file","var_cope_file")]),
    (inputnode,    level2model,   [(("cope", lambda x: len(x)), "num_copes")]),
    (level2model,  flameo,        [("design_mat","design_file"),
                                   ("design_con","t_con_file"),
                                   ("design_grp","cov_split_file")]),
    (flameo,       overlayflame,  [("zstats","stat_image")]),
    (overlayflame, sliceflame,    [("out_file", "in_file")]),
    (flameo,       outputnode,    [("stats_dir", "stats")]),
    (sliceflame,   report,        [("out_file", "zstat")]),
    ])

