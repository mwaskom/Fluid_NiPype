import nipype.algorithms.modelgen as modelgen
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe


# Model Workflow Definition
fsl_model = pe.Workflow(name="model")

inputnode = pe.Node(util.IdentityInterface(fields=["subject_info",
                                                   "TR", 
                                                   "units", 
                                                   "hpf_cutoff",
                                                   "HRF_bases", 
                                                   "contrasts",
                                                   "outlier_files",
                                                   "realignment_parameters",
                                                   "timeseries"]),
                    name="inputspec")

modelspec = pe.Node(modelgen.SpecifyModel(concatenate_runs=False),
                    name="modelspec")

level1design = pe.Node(fsl.Level1Design(model_serial_correlations="AR(1)"),
                       name="level1design")

modelgen = pe.MapNode(fsl.FEATModel(), name="modelgen",
                      overwrite = False,
                      iterfield = ["fsf_file"])

modelestimate = pe.MapNode(fsl.FILMGLS(smooth_autocorr=True,
                                       mask_size=5,
                                       threshold=1000),
                           name="modelestimate",
                           iterfield = ["design_file","in_file"])

threshres = pe.MapNode(fsl.ImageStats(op_string = "-r"),
                       name="threshresidual",
                       iterfield=["in_file"])
  
mni_brain = fsl.Info.standard_image("avg152T1.nii.gz")

overlayres = pe.MapNode(fsl.Overlay(auto_thresh_bg=True,
                                    background_image=mni_brain),
                        name="overlayresidual",
                        iterfield = ["stat_image", "stat_thresh"])

sliceres = pe.MapNode(fsl.Slicer(image_width=872),
                      name="sliceresidual",
                      iterfield=["in_file"])
sliceres.inputs.sample_axial = 2

contrastestimate = pe.MapNode(fsl.ContrastMgr(), name="contrastestimate",
                              iterfield = ["tcon_file","stats_dir"])

selectcontrast = pe.MapNode(util.Select(), 
                            name="selectcontrast", 
                            iterfield=["inlist"])

overlaystats = pe.MapNode(fsl.Overlay(stat_thresh=(2.3,10),
                                      auto_thresh_bg=True,
                                      show_negative_stats=True,
                                      background_image=mni_brain),
                          name="overlaystats",
                          iterfield = ["stat_image"])
                     
slicestats = pe.MapNode(fsl.Slicer(image_width=872),
                        name="slicestats",
                        iterfield=["in_file"])
slicestats.inputs.sample_axial = 2

outputnode = pe.Node(util.IdentityInterface(fields=["results",
                                                    "copes",
                                                    "varcopes",
                                                    "zstats"]),
                     name="outputspec")

report = pe.Node(util.IdentityInterface(fields=["design_image",
                                                "design_covariance",
                                                "residual",
                                                "zstat"]),
                 name="report")


def con_sort(files):
    files.sort()
    return files

fsl_model.connect([
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
    (modelestimate,     overlayres,        [("sigmasquareds","stat_image")]),
    (threshres,         overlayres,      
     [(("out_stat", lambda x: [tuple(i) for i in x]), "stat_thresh")]),
    (overlayres,        sliceres,          [("out_file", "in_file")]),
    (modelestimate,     contrastestimate,  [("results_dir","stats_dir")]),
    (selectcontrast,    overlaystats,      [("out","stat_image")]),
    (overlaystats,      slicestats,        [("out_file", "in_file")]),
    (modelestimate,     outputnode,        [("results_dir", "results")]),
    (contrastestimate,  outputnode,        [("copes", "copes"),
                                            ("varcopes", "varcopes"),
                                            ("zstats", "zstats")]),
    (modelgen,          report,            [("design_image", "design_image"),
                                            ("design_cov", "design_covariance")]),
    (sliceres,          report,            [("out_file", "residual")]),
    (slicestats,        report,            [("out_file", "zstat")]),
    ])

