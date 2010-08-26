import nipype.algorithms.modelgen as model
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe


fsl_modelfit = pe.Workflow(name='fsl_modelfit')

modelspec = pe.Node(interface=model.SpecifyModel(),  name="modelspec")
modelspec.inputs.concatenate_runs = False

level1design = pe.Node(interface=fsl.Level1Design(), name="level1design")

modelgen = pe.MapNode(interface=fsl.FEATModel(), name='modelgen',
                      iterfield = ['fsf_file'])

modelgen.overwrite = False

modelestimate = pe.MapNode(interface=fsl.FILMGLS(smooth_autocorr=True,
                                                 mask_size=5,
                                                 threshold=1000),
                           name='modelestimate',
                           iterfield = ['design_file','in_file'])

threshres = pe.MapNode(interface=fsl.ImageStats(op_string = '-r'),
                       name='threshresidual',
                       iterfield=['in_file'])

overlayres = pe.MapNode(interface=fsl.Overlay(auto_thresh_bg=True),
                        name='overlayresidual',
                        iterfield = ['stat_image', 'stat_thresh'])

sliceres = pe.MapNode(interface=fsl.Slicer(all_axial=True,
                                           image_width=750),
                      name='sliceresidual',
                      iterfield=['in_file'])

contrastestimate = pe.MapNode(interface=fsl.ContrastMgr(), name='contrastestimate',
                         iterfield = ['tcon_file','stats_dir'])

selectcontrast = pe.MapNode(interface=util.Select(), name='selectcontrast', iterfield=['inlist'])

overlaystats = pe.MapNode(interface=fsl.Overlay(stat_thresh=(2.5,10),
                                                auto_thresh_bg=True,
                                                show_negative_stats=True),
                         name='overlaystats',
                         iterfield = ['stat_image'])
                     
slicestats = pe.MapNode(interface=fsl.Slicer(all_axial=True,
                                             image_width=750),
                        name='slicestats',
                        iterfield=['in_file'])


fsl_modelfit.connect([
   (modelspec,level1design,[('session_info','session_info')]),
   (level1design,modelgen,[('fsf_files','fsf_file')]),
   (modelgen,modelestimate,[('design_file','design_file')]),
   (modelgen,contrastestimate,[('con_file','tcon_file')]),
   (modelestimate,threshres,[('sigmasquareds','in_file')]),
   (modelestimate,overlayres,[('sigmasquareds','stat_image')]),
   (threshres,overlayres,[(('out_stat', lambda x: [tuple(i) for i in x]), 'stat_thresh')]),
   (overlayres,sliceres,[('out_file', 'in_file')]),
   (modelestimate,contrastestimate,[('results_dir','stats_dir')]),
   (contrastestimate,selectcontrast,[('zstats','inlist')]),
   (selectcontrast, overlaystats, [('out','stat_image')]),
   (overlaystats, slicestats, [('out_file', 'in_file')]),
   ])

