from datetime import datetime

def subject_container(workflow, subjectsource, datasinknode, stripstring=None):
    
    if stripstring is None:
        outputs = subjectsource.outputs.get()
        stripstring = "_%s_"%outputs.keys()[0]
    workflow.connect([
        (subjectsource, datasinknode, 
            [("subject_id", "container"),
            (("subject_id", join_strings, stripstring), "strip_dir")]),
            ])

def join_strings(x, y):
    return "".join([y, x])


def substitute(origpath, subname):
    """Generate a list of substitution tuples."""
    import os
    from nipype.utils.filemanip import split_filename
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
        substitutes.append((os.path.basename(path), "".join([subname,ext])))
    return substitutes

def get_output_substitutions(workflow, outputnode, mergenode):
    """Substitute the output field name for the filename for all outputs from a node
    and send into a mergenode."""
    outputs = outputnode.outputs.get()
    for i, field in enumerate(outputs):
        workflow.connect(outputnode, (field, substitute, field), mergenode, "in%d"%(i+1))
    
def get_mapnode_substitutions(workflow, output_node, nruns):
    import networkx as nx
    from nipype.pipeline.engine import MapNode
    substitutions = []
    mapnodes = [e[0].name for e in nx.edges(workflow._graph) \
                    if e[1] is output_node and isinstance(e[0],MapNode)]
    for r in range(nruns):
        for node in mapnodes:
            substitutions.append(("_%s%d"%(node, r), "run_%d"%(r+1)))
    return substitutions

def set_substitutions(workflow, sinknode, mergenode, substitutions):

    workflow.connect(
    	mergenode, ("out", build_sub_list, substitutions), sinknode, "substitutions")

def build_sub_list(merged_subs, additional_subs):

    full_subs = merged_subs + additional_subs
    return full_subs

def connect_inputs(workflow, datagrabber, inputnode, makelist=[], listlength=None):
    """Connect the outputs of a Datagrabber to an inputnode.

    The names of the datagrabber outfields and the inputnode fields must match.

    The makelist parameter can be used to wrap certain inputs in a list.
    """
    inputs = inputnode.inputs.get()
    outputs = datagrabber.outputs.get()
    fields = [f for f in inputs if f in outputs]
    for field in fields:
        if field in makelist:
            workflow.connect(datagrabber, (field, make_list, listlength), inputnode, field)
        else:
            workflow.connect(datagrabber, field, inputnode, field)

def make_list(inputs, listlength):
    
    if listlength is not None and len(inputs) != listlength:
        listlength = [inputs for i in range(listlength)]
    elif not isinstance(inputs, list):
        inputs = [inputs]
    return inputs

def sink_outputs(workflow, outputnode, datasinknode, pathstr):
    
    if pathstr.endswith("."):
        pathstr = pathstr + "@"
    elif not pathstr.endswith(".@"):
        pathstr = pathstr + ".@"
    
    outputs = outputnode.outputs.get()
    for field in outputs:
        workflow.connect(outputnode, field, datasinknode, pathstr + field)

def archive_crashdumps(workflow):
    """Archive crashdumps by date to Nipype_Code directory"""
    import os
    datestamp = str(datetime.now())[:10]
    codepath = os.path.split(os.path.abspath(__file__))[0]
    crashdir = os.path.abspath("%s/crashdumps/%s" % (codepath, datestamp))
    if not os.path.isdir(crashdir):    
        os.makedirs(crashdir)
    workflow.config = dict(crashdump_dir=crashdir) 


def parse_par_file(parfile):
    onsets = []
    durations = []
    amplitudes = []
    for line in open(parfile):
        line = line.split()
        onsets.append(float(line[0]))
        durations.append(float(line[1]))
        amplitudes.append(float(line[2]))
    return onsets, durations, amplitudes

