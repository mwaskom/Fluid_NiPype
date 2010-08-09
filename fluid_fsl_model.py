import nipype.algorithms.modelgen as model
import nipype.interfaces.fsl as fsl
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

contrastestimate = pe.MapNode(interface=fsl.ContrastMgr(), name='contrastestimate',
                         iterfield = ['tcon_file','stats_dir'])

fsl_modelfit.connect([
   (modelspec,level1design,[('session_info','session_info')]),
   (level1design,modelgen,[('fsf_files','fsf_file')]),
   (modelgen,modelestimate,[('design_file','design_file')]),
   (modelgen,contrastestimate,[('con_file','tcon_file')]),
   (modelestimate,contrastestimate,[('results_dir','stats_dir')]),
   ])

