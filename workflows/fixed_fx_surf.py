import nipype.interfaces.fsl as fsl
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe
import numpy as np

# Inputs
inputnode = pe.Node(util.IdentityInterface(fields=["hemi","cope", "varcope"]),
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

# Just got a volume that (should be) all 1's to use as a mask for now
getmask = pe.Node(fs.Binarize(min=-np.inf), name="getmask")

# Run a fixed effects analysis in FLAMEO
flameo = pe.Node(fsl.FLAMEO(run_mode="fe"),
                 name="flameo")

# Take a picture of the lovely activations
snapshots = pe.Node(fs.SurfaceSnapshots(surface="cool",
                                        subject_id="fsaverage",
                                        overlay_range=(2.3,7),
                                        show_gray_curv=True,
                                        show_color_scale=True),
                    name="snapshots")

# Outputs
outputnode = pe.Node(util.IdentityInterface(fields=["stats"]), name="outputspec")

# Report
report = pe.Node(util.IdentityInterface(fields=["snapshots"]), name="report")

fixed_fx = pe.Workflow(name="fixed_fx")

fixed_fx.connect([
    (inputnode,    copemerge,     [("cope", "in_files")]),
    (inputnode,    varcopemerge,  [("varcope", "in_files")]),
    (copemerge,    flameo,        [("merged_file", "cope_file")]),
    (varcopemerge, flameo,        [("merged_file", "var_cope_file")]),
    (inputnode,    level2model,   [(("cope", lambda x: len(x)), "num_copes")]),
    (level2model,  flameo,        [("design_mat", "design_file"),
                                   ("design_con", "t_con_file"),
                                   ("design_grp", "cov_split_file")]),
    (inputnode,    getmask,       [(("cope", lambda x: x[0]), "in_file")]),
    (getmask,      flameo,        [("binary_file", "mask_file")]),
    (inputnode,    snapshots,     [("hemi", "hemi")]),
    (flameo,       snapshots,     [("zstats", "overlay")]),
    (flameo,       outputnode,    [("stats_dir", "stats")]),
    (snapshots,    report,        [("snapshots",  "snapshots")]),
    ])

