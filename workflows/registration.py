import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util


def get_registration_workflow(name="registration", volume=True, surface=True, surface_smooth=True):

    registration = pe.Workflow(name=name)

    # Define the inputs for the registation workflow
    infields = []
    if volume:
        infields.extend(["vol_source", "warpfield", "fsl_affine"])
    if surface:
        infields.extend(["surf_source", "subject_id", "tkreg_affine"])
    if surface_smooth:
        infields.extend(["smooth_fwhm"])
    inputnode = pe.Node(util.IdentityInterface(fields=infields),
                        name="inputspec")

    if volume:

        mni152 = fsl.Info.standard_image("avg152T1_brain.nii.gz")
        applywarp = pe.MapNode(fsl.ApplyWarp(ref_file=mni152),
                                iterfield=["in_file", "premat"],
                             name="applywarp")
        
        registration.connect([
            (inputnode, applywarp, [("vol_source", "in_file"),
                                    ("warpfield", "field_file"),
                                    ("fsl_affine", "premat")]),
                ])

    if surface:
        
        hemisource = pe.Node(util.IdentityInterface(fields=["hemi"]),
                             iterables=("hemi",["lh","rh"]),
                             name="hemisource")

        surfproject = pe.MapNode(fs.SampleToSurface(sampling_range=(0,1,.1),
                                                    sampling_units="frac",
                                                    cortex_mask=True),
                                 iterfield=["source_file", "reg_file"],
                                 name="surfproject")
        surfproject.inputs.sampling_method="average"

        surftransform = pe.MapNode(fs.SurfaceTransform(target_subject="fsaverage"),
                                   iterfield=["source_file"],
                                   name="surftransform")

        cvtnormsurf = pe.MapNode(fs.MRIConvert(out_type="niigz"),
                                 iterfield=["in_file"],
                                 name="convertnormsurf")
        
        registration.connect([
            (inputnode,    surfproject,    [("surf_source", "source_file"),
                                            ("subject_id", "subject_id"),
                                            ("tkreg_affine", "reg_file")]),
            (surfproject,  surftransform,  [("out_file", "source_file")]),
            (inputnode,    surftransform,  [("subject_id", "source_subject")]),
            (hemisource,   surftransform,  [("hemi", "hemi")]),
            ])

    if surface_smooth:

        smoothnormsurf = pe.MapNode(fs.SurfaceSmooth(subject_id="fsaverage",
                                                     reshape=True),
                                    iterfield=["in_file"],
                                    name="smoothnormsurf")

        smoothnatsurf = pe.MapNode(fs.SurfaceSmooth(),
                                   iterfield=["in_file"],
                                   name="smoothnativesurf")

        registration.connect([
            (surftransform,smoothnormsurf, [("out_file", "in_file")]),
            (hemisource,   smoothnormsurf, [("hemi", "hemi")]),
            (inputnode,    smoothnormsurf, [("smooth_fwhm", "fwhm")]),
            (surfproject,  smoothnatsurf,  [("out_file", "in_file")]),
            (hemisource,   smoothnatsurf,  [("hemi", "hemi")]),
            (inputnode,    smoothnatsurf,  [("subject_id", "subject_id"),
                                            ("smooth_fwhm", "fwhm")]),
            (smoothnormsurf, cvtnormsurf,  [("out_file", "in_file")]),
            ])
    elif surface:
        registration.connect(surftransform, "out_file", cvtnormsurf, "in_file")

    outfields = []
    if volume:
        outfields.append("warped_image")
    if surface:
        outfields.extend(["hemi_image", "hemi_image_fsaverage", "hemi"])
    outputnode = pe.Node(util.IdentityInterface(fields=["warped_timeseries",
                                                        "warped_example_func",
                                                        "warped_functional_mask",
                                                        "warped_mean_func",
                                                        "hemi_timeseries",
                                                        "hemi_timeseries_fsaverage",
                                                        "register",
                                                        "hemi"]),
                         name="outputspec")

    
    return registration
