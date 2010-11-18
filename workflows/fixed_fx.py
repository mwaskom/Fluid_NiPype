import os
import nipype.interfaces.io as nio
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.fsl as fsl
import nipype.algorithms.modelgen as model
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe

"""
Set up fixed-effects workflow
-----------------------------
"""

"""
Define all of the nodes
"""

"""
Use :class:`nipype.interfaces.fsl.Merge` to merge the copes and
varcopes for each condition
"""

copemerge = pe.MapNode(interface=fsl.Merge(dimension='t'),
                       iterfield=["in_files"],
                       name="copemerge")

varcopemerge = pe.MapNode(interface=fsl.Merge(dimension='t'),
                          iterfield=["in_files"],
                          name="varcopemerge")

"""
Use :class:`nipype.interfaces.fsl.L2Model` to generate subject and condition
specific level 2 model design files
"""

level2model = pe.Node(interface=fsl.L2Model(),
                      name='l2model')

"""
Use :class:`nipype.interfaces.fsl.FLAMEO` to estimate a second level model
"""

brain_mask = fsl.Info.standard_image("avg152T1_brain.nii.gz")

flamemask = pe.Node(fs.Binarize(min=10,in_file=brain_mask),
                    name="flamemask")

flameo = pe.MapNode(interface=fsl.FLAMEO(run_mode='fe'),
                 iterfield=["cope_file","var_cope_file"],
                 name="flameo")

mni_brain = fsl.Info.standard_image("avg152T1.nii.gz")
overlayflame = pe.MapNode(interface=fsl.Overlay(stat_thresh=(2.5, 10),
                                                auto_thresh_bg=True,
                                                show_negative_stats=True,
                                                background_image=mni_brain),
                          iterfield=["stat_image"],
                          name='overlayflame')

sliceflame = pe.MapNode(interface=fsl.Slicer(all_axial=True,
                                             image_width=750),
                        iterfield=["in_file"],
                        name='sliceflame')

"""
Connect the nodes in the volume workflow
"""

vol_fixed_fx = pe.Workflow(name='vol_fixed_fx')

vol_fixed_fx.connect([(copemerge,flameo,[('merged_file','cope_file')]),
                      (varcopemerge,flameo,[('merged_file','var_cope_file')]),
                      (level2model,flameo, [('design_mat','design_file'),
                                            ('design_con','t_con_file'),
                                            ('design_grp','cov_split_file')]),
                      (flamemask,flameo, [("binary_file", "mask_file")]),
                      (flameo,overlayflame,[('zstats','stat_image')]),
                      (overlayflame,sliceflame,[('out_file', 'in_file')])
                  ])

"""
Clone the volume workflow for the surface flow
"""

surf_fixed_fx = vol_fixed_fx.clone("surf_fixed_fx")
