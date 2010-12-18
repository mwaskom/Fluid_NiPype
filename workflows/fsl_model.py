import os
import nipype.algorithms.modelgen as modelgen
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe


inputnode = pe.Node(util.IdentityInterface(fields=["subject_info",
                                                   "TR", 
                                                   "units", 
                                                   "hpf_cutoff",
                                                   "HRF_bases", 
                                                   "contrasts",
                                                   "outlier_files",
                                                   "overlay_background",
                                                   "realignment_parameters",
                                                   "timeseries"]),
                    name="inputspec")

# Generate Nipype-style model information
modelspec = pe.Node(modelgen.SpecifyModel(concatenate_runs=False),
                    name="modelspec")

# Generate FSL-style model information
level1design = pe.Node(fsl.Level1Design(model_serial_correlations="AR(1)"),
                       name="level1design")

# Use model information to create fsf files
modelgen = pe.MapNode(fsl.FEATModel(), name="modelgen",
                      overwrite = False,
                      iterfield = ["fsf_file"])

# Use film_gls to estimate the model
modelestimate = pe.MapNode(fsl.FILMGLS(smooth_autocorr=True,
                                       mask_size=5,
                                       threshold=1000),
                           name="modelestimate",
                           iterfield = ["design_file","in_file"])

# Get the robust threshold of the sigmasquareds volume
threshres = pe.MapNode(fsl.ImageStats(op_string = "-r"),
                       name="threshresidual",
                       iterfield=["in_file"])

# Overlay the sigmasquareds in color
overlayres = pe.MapNode(fsl.Overlay(auto_thresh_bg=True),
                        name="overlayresidual",
                        iterfield = ["stat_image", "stat_thresh","background_image"])

# Slice the sigmasquareds into a png
sliceres = pe.MapNode(fsl.Slicer(),
                      name="sliceresidual",
                      iterfield=["in_file"])

# Estimate contrasts from the parameter effect size images
contrastestimate = pe.MapNode(fsl.ContrastMgr(), name="contrastestimate",
                              iterfield = ["tcon_file","stats_dir"])

# This node will iterate over each contrast for reporting
# (The iterables must be set elsewhere, as this workflow is 
# agnostic to model information)
selectcontrast = pe.MapNode(util.Select(), 
                            name="selectcontrast", 
                            iterfield=["inlist"])

# Overlay the zstats
overlaystats = pe.MapNode(fsl.Overlay(stat_thresh=(2.3,10),
                                      auto_thresh_bg=True,
                                      show_negative_stats=True),
                          name="overlaystats",
                          iterfield = ["stat_image","background_image"])

# Slice the zstats for reporting
slicestats = pe.MapNode(fsl.Slicer(),
                        name="slicestats",
                        iterfield=["in_file"])

# Define the workflow outputs
outputnode = pe.Node(util.IdentityInterface(fields=["results",
                                                    "copes",
                                                    "varcopes",
                                                    "zstats"]),
                     name="outputspec")

# Define the reporting outputs
report = pe.Node(util.IdentityInterface(fields=["design_image",
                                                "design_covariance",
                                                "residual",
                                                "zstat"]),
                 name="report")


def con_sort(files):
    """Take a list, sort it, and return it."""
    files.sort()
    return files

def get_sampling_rate(bg_image):
    """Sample overlay images every 2 slices if in MNI space, otherwise show every slice."""
    if isinstance(bg_image, list):
        bg_image = bg_image[0]
    try:
        # This heurstic is not perfect, but will do for us
        if bg_image.startswith(os.environ["FSLDIR"]):
            return 2
    except KeyError:
        return 1
    return 1

def get_image_width(bg_image):
    """Set the image width of the slicer png based on what space the background image is in."""
    if isinstance(bg_image, list):
        bg_image = bg_image[0]
    try:
        # This heurstic is not perfect, but will do for us
        if bg_image.startswith(os.environ["FSLDIR"]):
            return 872
    except KeyError:
        return 750
    return 750


# Model Workflow Definition
volume_model = pe.Workflow(name="model")

# Connect up the model workflow
volume_model.connect([
    (inputnode,         modelspec,         [("subject_info", "subject_info"),
                                            ("TR", "time_repetition"),
                                            ("units", "input_units"),
                                            ("units", "output_units"),
                                            ("hpf_cutoff", "high_pass_filter_cutoff"),
                                            ("timeseries", "functional_runs"),
                                            ("outlier_files", "outlier_files"),
                                            ("realignment_parameters","realignment_parameters")]),
    (inputnode,         level1design,      [("contrasts", "contrasts"),
                                            ("TR", "interscan_interval"),
                                            ("HRF_bases", "bases")]),
    (inputnode,         modelestimate,     [("timeseries", "in_file")]),
    (modelspec,         level1design,      [("session_info","session_info")]),
    (level1design,      modelgen,          [("fsf_files","fsf_file")]),
    (modelgen,          modelestimate,     [("design_file","design_file")]),
    (modelgen,          contrastestimate,  [("con_file","tcon_file")]),
    (contrastestimate,  selectcontrast,    [(("zstats", con_sort),"inlist")]),
    (modelestimate,     threshres,         [("sigmasquareds","in_file")]),
    (inputnode,         overlayres,        [("overlay_background", "background_image")]),
    (modelestimate,     overlayres,        [("sigmasquareds","stat_image")]),
    (threshres,         overlayres,        [(("out_stat", lambda x: [tuple(i) for i in x]), "stat_thresh")]),
    (inputnode,         sliceres,          [(("overlay_background", get_sampling_rate), "sample_axial"),
                                            (("overlay_background", get_image_width), "image_width")]),
    (overlayres,        sliceres,          [("out_file", "in_file")]),
    (modelestimate,     contrastestimate,  [("results_dir","stats_dir")]),
    (selectcontrast,    overlaystats,      [("out","stat_image")]),
    (inputnode,         overlaystats,      [("overlay_background", "background_image")]),
    (overlaystats,      slicestats,        [("out_file", "in_file")]),
    (inputnode,         slicestats,        [(("overlay_background", get_sampling_rate), "sample_axial"),
                                            (("overlay_background", get_image_width), "image_width")]),
    (modelestimate,     outputnode,        [("results_dir", "results")]),
    (contrastestimate,  outputnode,        [("copes", "copes"),
                                            ("varcopes", "varcopes"),
                                            ("zstats", "zstats")]),
    (modelgen,          report,            [("design_image", "design_image"),
                                            ("design_cov", "design_covariance")]),
    (sliceres,          report,            [("out_file", "residual")]),
    (slicestats,        report,            [("out_file", "zstat")]),
    ])

