import os
from copy import deepcopy
from nipype.interfaces.base import Bunch

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


