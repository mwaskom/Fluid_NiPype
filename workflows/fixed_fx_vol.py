import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe
import nibabel as nib


def get_fixedfx_workflow(name="fixed_fx"):

    # Define the workflow
    fixed_fx = pe.Workflow(name=name)

    # Set up the inputs
    inputnode = pe.Node(util.IdentityInterface(fields=["cope", 
                                                       "varcope",
                                                       "dof_file"]),
                        name="inputspec")

    # Use FSL's MNI space brainmask
    brain_mask = fsl.Info.standard_image("MNI152_T1_2mm_brain_mask.nii.gz")

    # Concatenate the Cope for each run
    copemerge = pe.Node(fsl.Merge(dimension="t"),
                        name="copemerge")

    # Concatenate the Varcope for each run
    varcopemerge = pe.Node(fsl.Merge(dimension="t"),
                           name="varcopemerge")
    
    # Get an image of the DOFs
    getdof = pe.MapNode(fsl.ImageMaths(suffix="_dof",
                                       in_file2=brain_mask),
                        iterfield=["in_file", "op_string"],
                        name="getdof")
    
    dofmerge = pe.Node(fsl.Merge(dimension="t"),
                       name="dofmerge")

    # Set up a FLAMEO model
    level2model = pe.Node(fsl.L2Model(),
                          name="l2model")

    # Run a fixed effects analysis in FLAMEO
    flameo = pe.Node(fsl.FLAMEO(run_mode="fe", 
                                mask_file=brain_mask),
                     name="flameo")

    # Display on the FSL template brain
    mni_brain = fsl.Info.standard_image("avg152T1_brain.nii.gz")

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
    outputnode = pe.Node(util.IdentityInterface(fields=["stats",
                                                        "zstat"]),
                         name="outputspec")

    def get_dof_opstring(doffile):
        with open(doffile) as f:
            dof = f.read().split()
        return "-mul 0 -add %d -mas"%dof

    fixed_fx.connect([
        (inputnode,    copemerge,     [("cope", "in_files")]),
        (inputnode,    varcopemerge,  [("varcope", "in_files")]),
        (inputnode,    getdof,        [("cope", "in_file")]),
        (inputnode,    getdof,        [(("dof_file", get_dof_opstring), "op_string" )]),
        (getdof,       dofmerge,      [("out_file", "in_files")]),
        (copemerge,    flameo,        [("merged_file","cope_file")]),
        (varcopemerge, flameo,        [("merged_file","var_cope_file")]),
        (dofmerge,     flameo,        [("out_file", "dof_var_cope_file")]),
        (inputnode,    level2model,   [(("cope", lambda x: len(x)), "num_copes")]),
        (level2model,  flameo,        [("design_mat","design_file"),
                                       ("design_con","t_con_file"),
                                       ("design_grp","cov_split_file")]),
        (flameo,       overlayflame,  [("zstats","stat_image")]),
        (overlayflame, sliceflame,    [("out_file", "in_file")]),
        (flameo,       outputnode,    [("stats_dir", "stats")]),
        (sliceflame,   outputnode,    [("out_file", "zstat")]),
        ])

    return fixed_fx, inputnode, outputnode
