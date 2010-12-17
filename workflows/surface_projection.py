import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe

inputnode = pe.Node(util.IdentityInterface(fields=["subject_id",
                                                   "reg_matrix",
                                                   "cope",
                                                   "varcope",
                                                   "smooth_fwhm"]),
                    name="inputspec")

hemisource = pe.Node(util.IdentityInterface(fields=["hemi"]),
                     iterables=("hemi", ["lh","rh"]),
                     name="hemisource")

# Actually sample the values onto the surface
copesampler = pe.MapNode(fs.SampleToSurface(sampling_range=(0,1,.1),
                                           sampling_units="frac",
                                           out_type="niigz",
                                           cortex_mask=True),
                        iterfield=["source_file","reg_file"],
                        name="copesampler")
copesampler.inputs.sampling_method = "average"

varcopesampler = pe.MapNode(fs.SampleToSurface(sampling_range=(0,1,.1),
                                              sampling_units="frac",
                                              out_type="niigz",
                                              cortex_mask=True),
                           iterfield=["source_file","reg_file"],
                           name="varcopesampler")
varcopesampler.inputs.sampling_method = "average"


# Transform the copes and varcopes to the fsaverage surface
copexfm = pe.MapNode(fs.SurfaceTransform(target_subject="fsaverage"),
                     iterfield=["source_file"],
                     name="copetransform")

varcopexfm = pe.MapNode(fs.SurfaceTransform(target_subject="fsaverage"),
                        iterfield=["source_file"],
                        name="varcopetransform")
# Smooth on the surface
copesmoother = pe.MapNode(fs.SurfaceSmooth(subject_id="fsaverage",
                                           reshape=True),
                          iterfield=["in_file"],
                          name="copesmoother")

varcopesmoother = pe.MapNode(fs.SurfaceSmooth(subject_id="fsaverage",
                                              reshape=True),
                             iterfield=["in_file"],
                             name="varcopesmoother")

# Define the outputs
outputnode = pe.Node(util.IdentityInterface(fields=["cope", "varcope"]),
                     name="outputspec")

surfproj = pe.Workflow(name="surfprojection")

surfproj.connect([
    (inputnode,       copesampler,     [("subject_id", "subject_id"),
                                        ("reg_matrix", "reg_file"),
                                        ("cope", "source_file")]),
    (inputnode,       varcopesampler,  [("subject_id", "subject_id"),
                                        ("reg_matrix", "reg_file"),
                                        ("varcope", "source_file")]),
    (hemisource,      copesampler,     [("hemi", "hemi")]),
    (hemisource,      varcopesampler,  [("hemi", "hemi")]),
    (hemisource,      copexfm,         [("hemi", "hemi")]),
    (inputnode,       copexfm,         [("subject_id", "source_subject")]),
    (copesampler,     copexfm,         [("out_file", "source_file")]),
    (hemisource,      varcopexfm,      [("hemi", "hemi")]),
    (inputnode,       varcopexfm,      [("subject_id", "source_subject")]),
    (varcopesampler,  varcopexfm,      [("out_file", "source_file")]),
    (inputnode,       copesmoother,    [("smooth_fwhm", "fwhm")]),
    (hemisource,      copesmoother,    [("hemi", "hemi")]),
    (copexfm,         copesmoother,    [("out_file", "in_file")]),
    (inputnode,       varcopesmoother, [("smooth_fwhm", "fwhm")]),
    (hemisource,      varcopesmoother, [("hemi", "hemi")]),
    (varcopexfm,      varcopesmoother, [("out_file", "in_file")]),
    (copesmoother,    outputnode,      [("out_file", "cope")]),
    (varcopesmoother, outputnode,      [("out_file", "varcope")])
    ])


