
import os
import nipype.interfaces.io as nio
import nipype.interfaces.fsl as fsl
import nipype.algorithms.modelgen as model
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe


modelfit = pe.Workflow(name='modelfit')

"""
Use :class:`nipype.algorithms.modelgen.SpecifyModel` to generate design information.
"""

modelspec = pe.Node(interface=model.SpecifyModel(),  name="modelspec")
modelspec.inputs.concatenate_runs = False

"""
Use :class:`nipype.interfaces.fsl.Level1Design` to generate a run specific fsf
file for analysis
"""

level1design = pe.Node(interface=fsl.Level1Design(), name="level1design")

"""
Use :class:`nipype.interfaces.fsl.FEATModel` to generate a run specific mat
file for use by FILMGLS
"""

modelgen = pe.MapNode(interface=fsl.FEATModel(), name='modelgen',
                      iterfield = ['fsf_file'])

"""
Set the model generation to run everytime. Since the fsf file, which is the
input to modelgen only references the ev files, modelgen will not run if the ev
file contents are changed but the fsf file is untouched.
"""

modelgen.overwrite = True

"""
Use :class:`nipype.interfaces.fsl.FILMGLS` to estimate a model specified by a
mat file and a functional run
"""

modelestimate = pe.MapNode(interface=fsl.FILMGLS(smooth_autocorr=True,
                                                 mask_size=5,
                                                 threshold=1000),
                           name='modelestimate',
                           iterfield = ['design_file','in_file'])

"""
Use :class:`nipype.interfaces.fsl.ContrastMgr` to generate contrast estimates
"""

conestimate = pe.MapNode(interface=fsl.ContrastMgr(), name='conestimate',
                         iterfield = ['tcon_file','stats_dir'])

modelfit.connect([
   (modelspec,level1design,[('session_info','session_info')]),
   (level1design,modelgen,[('fsf_files','fsf_file')]),
   (modelgen,modelestimate,[('design_file','design_file')]),
   (modelgen,conestimate,[('con_file','tcon_file')]),
   (modelestimate,conestimate,[('results_dir','stats_dir')]),
   ])

