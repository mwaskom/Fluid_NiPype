import nipype.algorithms.modelgen as model
import nipype.interfaces.spm as spm
import nipype.pipeline.engine as pe


spm_modelfit = pe.Workflow(name= "spm_modelfit")

modelspec = pe.Node(interface=model.SpecifyModel(),  name="modelspec")
modelspec.inputs.concatenate_runs = False

level1design = pe.Node(interface=spm.Level1Design(), name= "level1design")

modelestimate = pe.Node(interface=spm.EstimateModel(estimation_method = {"Classical" : 1}),
                                                    name="modelestimate")

contrastestimate = pe.Node(interface = spm.EstimateContrast(), name="contrastestimate")


spm_modelfit.connect([(modelspec, level1design, [("session_info", "session_info")]),
                      (level1design, modelestimate, [("spm_mat_file", "spm_mat_file")]),
                      (modelestimate, contrastestimate,[('spm_mat_file','spm_mat_file'),
                                                        ('beta_images','beta_images'),
                                                        ('residual_image','residual_image')]),
                      ])
