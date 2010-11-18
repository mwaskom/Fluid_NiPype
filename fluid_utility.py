import os
from copy import deepcopy
from datetime import datetime
from nipype.interfaces.base import Bunch
from nipype.utils.filemanip import split_filename

def subject_container(workflow, subjectsource, datasinknode, stripstring=None):
    
    if stripstring is None:
        outputs = subjectsource.outputs.get()
        stripstring = "_%s_"%outputs.keys()[0]
    workflow.connect([
        (subjectsource, datasinknode, 
            [("subject_id", "container"),
            (("subject_id", lambda x: "".join([stripstring,x])), "strip_dir")]),
            ])

def get_substitutions(workflow, outputnode, mergenode):
    """Substitute the output field name for the filename for all outputs from a node
    and send into a mergenode."""
    outputs = outputnode.outputs.get()
    for i, field in enumerate(outputs):
        workflow.connect(outputnode, (field, substitute, field), mergenode, "in%d"%(i+1))
    
def substitute(origpath, subname):
    """Generate a list of substitution tuples."""
    if not isinstance(origpath, list):
        origpath = [origpath]
    substitutes = []
    for path in origpath:
        ext = split_filename(path)[2]
        # Text files (from ART, etc.) tend to have an image 
        # filename hanging around, so this solution is a bit
        # messier in code but gets us better filenames.
        if ext.startswith(".nii.gz") and not ext == ".nii.gz":
            ext = ext[7:]
        elif ext.endswith(".txt"):
            ext = ".txt"
        # .mincost is a dumb extension
        elif ext.endswith(".mincost"):
            ext = ".dat"
        substitutes.append((os.path.basename(path), subname + ext))
    return substitutes

def connect_inputs(workflow, datagrabber, inputnode):
    """Connect the outputs of a Datagrabber to an inputnode.

    The names of the datagrabber outfields and the inputnode fields must match.
    """
    inputs = inputnode.inputs.get()
    outputs = datagrabber.outputs.get()
    fields = [f for f in inputs if f in outputs]
    for field in fields:
        workflow.connect(datagrabber, field, inputnode, field)

def sink_outputs(workflow, outputnode, datasinknode, pathstr):
    
    if pathstr.endswith("."):
        pathstr = pathstr + "@"
    elif not pathstr.endswith(".@"):
        pathstr = pathstr + ".@"
    
    outputs = outputnode.outputs.get()
    for field in outputs:
        workflow.connect(outputnode, field, datasinknode, pathstr + field)

def archive_crashdumps(workflow):
    """Archive crashdumps by date to NiPype_Code directory"""
    datestamp = str(datetime.now())[:10]
    codepath = os.path.split(os.path.abspath(__file__))[0]
    crashdir = os.path.abspath("%s/crashdumps/%s" % (codepath, datestamp))
    if not os.path.isdir(crashdir):    
        os.makedirs(crashdir)
    workflow.config = dict(crashdump_dir=crashdir) 


def parse_par_file(parfile):
    onsets = []
    durations = []
    for line in open(parfile):
        line = line.split()
        onsets.append(float(line[0]))
        durations.append(float(line[1]))
    return onsets, durations

