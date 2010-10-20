import os
from copy import deepcopy
from datetime import datetime
from nipype.interfaces.base import Bunch
from nipype.utils.filemanip import split_filename

def subject_container(workflow, subjectsource, datasinknode, stripstring="_subject_id_"):
    
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
        elif ext.endswith(".mincost"):
            ext = ".dat"
        substitutes.append((os.path.basename(path), subname + ext))
    return substitutes

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


def sort_copes(files):
    numelements = len(files[0])
    outfiles = []
    for i in range(numelements):
        outfiles.insert(i,[])
        for j, elements in enumerate(files):
            outfiles[i].append(elements[i])
    return outfiles

def parse_par_file(base_dir, subject_id, run_number, name, template):
    parfile = os.path.join(base_dir, template % (subject_id, run_number, name, subject_id[-4:]))
    fid = open(parfile,"r")
    onsets = []
    durations = []
    for line in fid:
        line = line.split()
        onsets.append(float(line[0]))
        durations.append(float(line[1]))
    return onsets, durations

def subjectinfo(subject_id, exp):
    print "Subject ID: %s\n"%str(subject_id)
    output = []
    names = exp.names
    for r in range(exp.nruns):
        onsets = []
        durations = []
        for name in names:
            o, d = parse_par_file(exp.data_dir, subject_id, r+1, name, exp.parfile_template)
            onsets.append(o)
            durations.append(d)
        output.insert(r,
                      Bunch(conditions=names,
                            onsets=deepcopy(onsets),
                            durations=deepcopy(durations),
                            amplitudes=None,
                            tmod=None,
                            pmod=None,
                            regressor_names=None,
                            regressors=None))
    return output

