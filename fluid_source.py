"""
    Info and data source nodes for fluid project.
"""    

import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util

data_dir = "/mindhive/gablab/fluid/Data"

infosource = pe.Node(interface=util.IdentityInterface(fields=["subject_id"]),
                     name="infosource")


datasource = pe.Node(interface=nio.DataGrabber(infields=["subject_id"],
                                               outfields=["func", "target", "struct"],
                                               base_directory = data_dir),
                     name="datasource")

