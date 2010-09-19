import nipype.algorithms.modelgen as model
import nipype.algorithms.rapidart as art
import nipype.interfaces.fsl as fsl
import nipype.interfaces.spm as spm
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe


"""
First level model workflow in the volume
"""

modelspec = pe.Node(interface=model.SpecifyModel(concatenate_runs=False),
                    name="modelspec")

level1design = pe.Node(interface=fsl.Level1Design(), name="level1design")

modelgen = pe.MapNode(interface=fsl.FEATModel(), name="modelgen",
                      iterfield = ["fsf_file"])

modelgen.overwrite = False

modelestimate = pe.MapNode(interface=fsl.FILMGLS(smooth_autocorr=True,
                                                 mask_size=5,
                                                 threshold=1000),
                           name="modelestimate",
                           iterfield = ["design_file","in_file"])

threshres = pe.MapNode(interface=fsl.ImageStats(op_string = "-r"),
                       name="threshresidual",
                       iterfield=["in_file"])
  
overlayres = pe.MapNode(interface=fsl.Overlay(auto_thresh_bg=True),
                        name="overlayresidual",
                        iterfield = ["stat_image", "stat_thresh"])

sliceres = pe.MapNode(interface=fsl.Slicer(all_axial=True,
                                           image_width=750),
                      name="sliceresidual",
                      iterfield=["in_file"])

contrastestimate = pe.MapNode(interface=fsl.ContrastMgr(), name="contrastestimate",
                              iterfield = ["tcon_file","stats_dir"])

mergecontrast = pe.MapNode(interface=util.Merge(3,axis="hstack"),
                           iterfield = ["in1","in2","in3"],
                           name = "mergecontrast")

selectcontrast = pe.MapNode(interface=util.Select(), 
                            name="selectcontrast", 
                            iterfield=["inlist"])

xfmcopes = pe.MapNode(interface=fs.ApplyVolTransform(fs_target=True,
                                                     no_resample=True),
                      iterfield = ["source_file","reg_file"],
                      name = "xfmcopes")

xfmvarcopes = pe.MapNode(interface=fs.ApplyVolTransform(fs_target=True,
                                                        no_resample=True),
                         iterfield = ["source_file","reg_file"],
                         name = "xfmvarcopes")

overlaystats = pe.MapNode(interface=fsl.Overlay(stat_thresh=(2.5,10),
                                                auto_thresh_bg=True,
                                                show_negative_stats=True),
                          name="overlaystats",
                          iterfield = ["stat_image"])
                     
slicestats = pe.MapNode(interface=fsl.Slicer(all_axial=True,
                                            image_width=750),
                        name="slicestats",
                        iterfield=["in_file"])

"""
Connect all of the nodes in a volume workflow
"""

fsl_vol_model = pe.Workflow(name="fsl_vol_model")

def con_sort(files):
    files.sort()
    return files

fsl_vol_model.connect([
   (modelspec,level1design,[("session_info","session_info")]),
   (level1design,modelgen,[("fsf_files","fsf_file")]),
   (modelgen,modelestimate,[("design_file","design_file")]),
   (modelgen,contrastestimate,[("con_file","tcon_file")]),
   (contrastestimate,mergecontrast,[(("copes", con_sort), "in1"),
                                    (("varcopes", con_sort), "in2"),
                                    (("zstats", con_sort), "in3")]),
   (mergecontrast,selectcontrast,[("out","inlist")]),
   (selectcontrast,xfmcopes,[(("out",lambda x: [l[0] for l in x]),"source_file")]),
   (selectcontrast,xfmvarcopes,[(("out",lambda x: [l[1] for l in x]),"source_file")]),
   (modelestimate,threshres,[("sigmasquareds","in_file")]),
   (modelestimate,overlayres,[("sigmasquareds","stat_image")]),
   (threshres,overlayres,[(("out_stat", lambda x: [tuple(i) for i in x]), "stat_thresh")]),
   (overlayres,sliceres,[("out_file", "in_file")]),
   (modelestimate,contrastestimate,[("results_dir","stats_dir")]),
   (selectcontrast,overlaystats,[(("out",lambda x: [l[2] for l in x]),"stat_image")]),
   (overlaystats, slicestats, [("out_file", "in_file")]),
   ])

"""
Clone the volume modelfitting workflow for the surface stream
"""

fsl_surf_model = fsl_vol_model.clone("fsl_surf_model")
