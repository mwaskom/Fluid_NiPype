import os
from nipype.interfaces.base import TraitedSpec, File, traits
from nipype.interfaces.fsl.base import FSLCommand, FSLCommandInputSpec
from nipype.utils.filemanip import fname_presuffix


class CheckRegInput(FSLCommandInputSpec):
    
    in_file = File(exists=True, argstr="%s", position=1)
    out_file = File(genfile=True, argstr="%s", position=2)

class CheckRegOutput(TraitedSpec):

    out_file = File(exists=True)

class CheckReg(FSLCommand):

    _cmd = "check_mni_reg"
    input_spec = CheckRegInput
    output_spec = CheckRegOutput

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["out_file"] = fname_presuffix(self.inputs.in_file,
                                              suffix="_to_mni.png",
                                              use_ext=False,
                                              newpath=os.getcwd())
        return outputs

    def _gen_filename(self, name):
        if name == "out_file":
            return self._list_outputs()[name]
        return None

class TimeSeriesMovieInput(FSLCommandInputSpec):

    in_file = File(exists=True,argstr="-ts %s")
    ref_type = traits.String(argstr="-ref %s")
    plot_file = File(exists=True,argstr="-plot %s")
    norm_plot = traits.Bool(argstr="-normplot")
    art_min = traits.Float(argstr="-min %.3f")
    art_max = traits.Float(argstr="-max %.3f")
    out_file = traits.File(genfile=True, argstr="-out %s")

class TimeSeriesMovieOutput(TraitedSpec):

    out_file = File(exists=True)

class TimeSeriesMovie(FSLCommand):

    _cmd = "ts_movie"
    input_spec = TimeSeriesMovieInput
    output_spec = TimeSeriesMovieOutput

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["out_file"] = fname_presuffix(self.inputs.in_file,
                                              suffix=".gif",
                                              use_ext=False,
                                              newpath=os.getcwd())
        return outputs

    def _gen_filename(self, name):
        if name == "out_file":
            return self._list_outputs()[name]
        return None


