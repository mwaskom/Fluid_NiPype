#! /usr/bin/env python

import os
import sys
import shutil
from glob import glob
from os.path import join as pjoin
import numpy as np
from scipy import isnan
from scipy.io import loadmat
import nibabel as nib
from htmltools import HTMLReport


DATA_DIR = "/mindhive/gablab/fluid/Data"
ANALYSIS_DIR = "/mindhive/gablab/fluid/Analysis/Nipype"
REPORT_DIR = "/mindhive/gablab/fluid/Analysis/Report"
STAGES = ["behavioral", "timeseries", "preproc", "registration" ,"model", "stats", "fixed_effects"]
FUNC_RUNS = dict(iq=1,nback=4,mot_block=1,mot_jitter=1,resting=1)
PARADIGM_MAP = dict(iq="IQ",nback="NBack",mot_block="MOT_Block",mot_jitter="MOT_Jitter",resting="Resting",mot_jitter_speed="MOT_Jitter")

def main():

    if len(sys.argv) == 1:
        sys.exit("USAGE: build_report.py [subj ...]")
   
    subjects = [s for s in sys.argv[1:] if os.path.exists(os.path.join(DATA_DIR, s))]

    for subj in subjects:
    
        topics = ["recon-all","fnirt"]
        for paradigm in FUNC_RUNS:
            if os.path.exists(os.path.join(ANALYSIS_DIR, paradigm, subj)):
                topics.append(paradigm)

        for topic in topics:
            if topic in FUNC_RUNS:
                for stage in STAGES:
                    write_func_report(subj, topic, stage)

            write_topic_home(subj, topic)
            
        write_subject_home(subj)

    write_homepage()

def write_subject_home(subj):
    
    html = HTMLReport(os.path.join(DATA_DIR, subj, "report", "report.html"))

    html.write_title("Subject %s Report"%subj)
    write_back_link(html, "home")

    html.write_section_head("Reports")

    topics = []
    for topic in ["recon-all", "fnirt"] + FUNC_RUNS.keys():
        if os.path.exists(os.path.join(DATA_DIR, subj, "report", "%s_report.html"%topic)):
            topics.append(topic)
    html.write_link_table("%s_report.html", "%s", topics, 10)

def write_back_link(html, targ, subj=None):

    if targ == "home":
        html.write_link(
            "Back to homepage",os.path.join(REPORT_DIR, "fluid_report.html"), "head")
        html.newline()
    elif targ == "subject":
        html.write_link(
            "Back to subject page",os.path.join(DATA_DIR, subj, "report", "report.html"),"head")
        html.newline()

def write_topic_home(subj, topic):
    
    html = HTMLReport(os.path.join(DATA_DIR, subj, "report", "%s_report.html"%topic))
    write_topic_head(html, subj, topic)
    if topic == "recon-all":
        write_recon_report(html, subj)
    elif topic == "fnirt":
        write_fnirt_report(html, subj)
    elif topic == "dti_vbm":
        write_dti_report(html, subj)

def write_topic_head(html, subj, topic):
    
    html.write_title("Subject: %s<br>%s report"%(subj, topic))
    write_back_link(html, "subject", subj)
    
    if topic in FUNC_RUNS:
        write_func_links(html, subj, topic)

def write_func_report(subj, paradigm, stage):

    html = HTMLReport(os.path.join(DATA_DIR, subj, "report", "%s_%s.html"%(paradigm, stage)))
    write_topic_head(html, subj, paradigm)

    if stage == "behavioral":
        write_bhvl_report(html, subj, paradigm)
    elif stage == "timeseries":
        write_timeseries_report(html, subj, paradigm)
    elif stage == "preproc":
        write_preproc_report(html, subj, paradigm)
    elif stage == "registration":
        write_registration_report(html, subj, paradigm)
    if stage != "resting":
        if stage == "model":
            write_model_report(html, subj, paradigm)
        elif stage == "stats":
            write_stats_report(html, subj, paradigm)
        elif stage == "fixed_effects":
	    write_ffx_report(html, subj, paradigm)

def write_recon_report(html, subj):

    status_log = os.path.join(DATA_DIR, subj, "scripts","recon-all-status.log")
    if os.path.exists(status_log):
        status = [l for l in open(status_log)][-1]
        html.write_text("Recon status: %s"%status)
        version = open(os.path.join(DATA_DIR, subj, "scripts","build-stamp.txt")).read()
        html.newline()
        html.write_text("Freesurfer version: %s"%version)
        recon_log = os.path.join(DATA_DIR, subj, "scripts", "recon-all.log")
        recon_text = os.path.join(DATA_DIR, subj, "scripts", "recon_log.txt")
	samba_recon_text = os.path.join(DATA_DIR, subj, "scripts", "recon_log.txt")
        shutil.copy(recon_log, recon_text)
        html.write_link("Recon-all Log", samba_recon_text)
        html.write_text("Surface Placement",bold=True)
        html.newline()
        html.write_image(os.path.join(DATA_DIR, subj, "gifs", "surf_movie.gif"))
        html.newline()
        html.write_text("Subcortical Segmentation",bold=True)
        html.newline()
        html.write_image(os.path.join(DATA_DIR, subj, "gifs", "aseg_movie.gif"))
        html.newline()
        html.write_text("White Matter Segmentation",bold=True)
        html.newline()
        html.write_image(os.path.join(DATA_DIR, subj, "gifs", "wmparc_movie.gif"))
    else:
        html.write_text("No recon-all status log exists")

def write_fnirt_report(html, subj):
    
    html.write_text("Final FNIRT Registration",bold=True)
    html.newline()
    html.write_image(os.path.join(DATA_DIR, subj, "normalization", "brain_warp_to_mni.png"))

def write_dti_report(html, subj):
    
    html.write_text("FA Map Normalization", bold=True)
    html.newline()
    html.write_image(os.path.join(
        ANALYSIS_DIR,"dti_vbm","report",subj,"preproc","reg_slices.png"))
    html.newline()
    html.write_image(os.path.join(
        ANALYSIS_DIR,"dti_vbm","report",subj,"preproc","image_slices.png"))

def write_func_links(html, subj, paradigm):

    if paradigm == "resting":
        func_stages = STAGES[1:3]
    else:
        func_stages = STAGES
    html.write_link_table(paradigm +"_%s.html", "%s", func_stages, 10)

def write_bhvl_report(html, subj, paradigm):

    srcdir = os.path.join(DATA_DIR, subj, "fmri_bhvl")
    if paradigm == "nback":
        write_nback_bhvl(html, subj, srcdir)
    elif paradigm == "iq":
        write_iq_bhvl(html, subj, srcdir)


def write_iq_bhvl(html, subj, srcdir):

    html.newline()
    try:
        ptbfile = glob(os.path.join(srcdir, "*_BLOCK_CATTELL.mat"))[0]
        stE = loadmat(ptbfile,squeeze_me=True,struct_as_record=False)["stE"]
    except IndexError:
        html.write_text("Could not open behavioral data file")
        html.newline(2)
        return
    
    diffdict = {1:"easy",0:"hard"}
    resps = dict(easy=0,hard=0)
    accs = dict(easy=[],hard=[])
    rt = dict(easy=[],hard=[])

    for t, acc in enumerate(stE.Trials.Accuracy):
        diff = diffdict[stE.Trials.bLowG[t]]
        if not isnan(acc):
            resps[diff] += 1
            accs[diff].append(acc)
            rt[diff].append(stE.Trials.RT[t])
    
    for diff in ["hard","easy"]:
        html.write_text("%s Stimuli"%diff.capitalize(), bold=True)
        html.newline()
        html.write_text("Responses: %d"%resps[diff])
        html.newline()
        html.write_text("Accuracy: %.3f"%np.mean(accs[diff]))
        html.newline()
        html.write_text("Mean RT: %.3f"%np.mean(rt[diff]))
        html.newline(2)


def write_nback_bhvl(html, subj, srcdir):
       
    html.newline()
    sessaccall = {}
    sessacchit = {}
    sessaccnon = {}
    for level in range(4):
        sessaccall[level] = []
        sessacchit[level] = []
        sessaccnon[level] = []
    for r in range(1,5):
        try:
            ptbfile = glob(os.path.join(srcdir, "*_run%d_DUALNBACK.mat"%r))[0]
        except IndexError:
            html.write_text("Could not open run %d data file"%r)
            html.newline(2)
            continue
        stE = loadmat(ptbfile,struct_as_record=False,squeeze_me=True)["stE"]
        html.write_text("Run %d Accuracy:"%r,bold=True)
        accstring = ""
        allaccstring = "<br>All Trials<br>&nbsp;&nbsp;"
        hitaccstring = "<br>Hit Trials<br>&nbsp;&nbsp;"
        nonaccstring = "<br>Null Trials<br>&nbsp;&nbsp;"
        for level in range(4):
            allacc = np.mean([np.mean(stE.Block[i].Trials.iAccuracy) for i in 
                            range(len(stE.Block)) if stE.Block[i].iNumBack == level])
            allaccstring += "%d-Back: %.3f"%(level, allacc)
            allaccstring += "".join("&nbsp;" for i in range(5))
            sessaccall[level].append(allacc)
            hitacc = np.mean([np.mean(stE.Block[i].Trials.iAccuracy[np.where(stE.Block[i].Trials.sCorrectResponse)]) 
                           for i in range(len(stE.Block)) if stE.Block[i].iNumBack == level])
            hitaccstring += "%d-Back: %.3f"%(level, hitacc)
            hitaccstring += "".join("&nbsp;" for i in range(5))
            sessacchit[level].append(hitacc)
            # If this makes sense next time I look at it, I get a cookie
            nonacc = np.mean([np.mean(stE.Block[j].Trials.iAccuracy[
                 np.array([i for i, val in enumerate(stE.Block[j].Trials.sCorrectResponse) 
                 if val == np.array([])],dtype=np.int)])
                 for j in range(len(stE.Block)) if stE.Block[j].iNumBack == level])
            nonaccstring += "%d-Back: %.3f"%(level, nonacc)
            nonaccstring += "".join("&nbsp;" for i in range(5))
            sessaccnon[level].append(nonacc)

        html.write_text(allaccstring)
        html.write_text(hitaccstring)
        html.write_text(nonaccstring)
        html.newline(2)
    html.write_text("Session accuracy:",bold=True)
    html.newline()
    allaccstring = "All Trials<br>&nbsp;&nbsp;"
    hitaccstring = "<br>Hit Trials<br>&nbsp;&nbsp;"
    nonaccstring = "<br> Null Trials<br>&nbsp;&nbsp;"
    for level in range(4):
        allaccstring += "%d-Back: %.3f"%(level, np.mean(sessaccall[level]))
        allaccstring += "".join("&nbsp;" for i in range(5))
        hitaccstring += "%d-Back: %.3f"%(level, np.mean(sessacchit[level]))
        hitaccstring += "".join("&nbsp;" for i in range(5))
        nonaccstring += "%d-Back: %.3f"%(level, np.mean(sessaccnon[level]))
        nonaccstring += "".join("&nbsp;" for i in range(5))
    html.write_text(allaccstring)
    html.write_text(hitaccstring)
    html.write_text(nonaccstring)
    
def write_timeseries_report(html, subj, paradigm):

    for r in range(1, FUNC_RUNS[paradigm] + 1):
        
        html.newline()
        html.write_section_head("Run %d"%r)
        html.newline()
        
        if paradigm != "resting":
            tsfile = os.path.join(DATA_DIR, subj, "bold", "%s_run%d.nii.gz"%(PARADIGM_MAP[paradigm], r))
        else:
            tsfile = os.path.join(DATA_DIR, subj, "bold", "%s.nii.gz"%(PARADIGM_MAP[paradigm]))

        html.write_text("Source file: %s"%tsfile)
        html.newline()
        try:
            html.write_text("Original path: %s"%os.path.realpath(tsfile))
            html.newline()
        except OSError:
            pass
        try:
            img = nib.load(tsfile)
            html.write_text("Image Dimensions: %dx%dx%d"%img.get_shape()[:3])
            html.newline()
            html.write_text("Timepoints: %d"%img.get_shape()[-1])
            html.newline()
        except IOError:
            html.write_text("Source timeseries could not be read by Nibabel")

        moviefile = os.path.join(ANALYSIS_DIR, paradigm, subj, "preproc", "run_%d"%r, "timeseries_movie.gif")
        html.write_image(moviefile)

def write_preproc_report(html, subj, paradigm):

    for r in range(1, FUNC_RUNS[paradigm] + 1):
    
        srcdir = os.path.join(ANALYSIS_DIR, paradigm, subj, "preproc", "run_%d"%r)

        html.newline()
        html.write_section_head("Run %d"%r)
        html.newline()

        html.write_text("Motion Correction Target",bold=True)
        html.write_image(os.path.join(srcdir, "example_func_slices.png"))
        html.newline()
        
        html.write_text("Mean Functional Image",bold=True)
        html.write_image(os.path.join(srcdir, "mean_func_slices.png"))
        html.newline()

        html.write_text("Mean Intensity Plot",bold=True)
        html.write_image(os.path.join(srcdir, "intensity_plot.png"))
        html.newline()

        try:
            outfile = os.path.join(ANALYSIS_DIR,paradigm,subj,"preproc","run_%d"%r,"outlier_volumes.txt")
            nout = len(open(outfile).read().strip().split())
            html.write_text("Total Outlier Volumes:&nbsp;",bold=True)
            html.write_text(nout)
        except IOError:
            html.write_text("Could not open %s for reading"%outfile)
        html.newline(2)

        try:
            maxfile =  pjoin(ANALYSIS_DIR,paradigm,subj,"preproc","run_%d"%r,"max_motion.txt")
            info = [l.strip() for l in open(maxfile).readlines()]
            html.write_text("Maximum RMS Motion:", bold=True)
            html.newline()
            html.write_text("Absolute -- %s mm"%info[1])
            html.newline()
            html.write_text("Relative -- %s mm"%info[3])
        except IOError:
            html.write_text("Could not open %s for reading"%maxfile)
        html.newline(2)

        html.write_text("Motion Plots",bold=True)

        for plot_type in ["rotation", "translation", "displacement"]:
            image = os.path.join(srcdir, "%s_plot.png"%plot_type) 
            html.write_image(image)

def write_registration_report(html, subj, paradigm):

    for r in range(1, FUNC_RUNS[paradigm] + 1):
        
        html.write_section_head("Run %s"%r)

        srcdir = os.path.join(ANALYSIS_DIR, paradigm, subj, "preproc", "run_%d"%r)
        
        costfile = os.path.join(srcdir, "func2anat.mincost")
        try:
            cost = open(costfile).read().split()[0]
            html.write_text("Final boundary-based registration cost value: %s"%cost)
        except IOError:
            html.write_text("Could not open %s for reading"%costfile)
        html.newline(2)

        html.write_text("Example func to native anatomical")
        image = os.path.join(srcdir, "func2anat_slices.png")
        html.write_image(image)
        html.newline()

def write_model_report(html, subj, paradigm):

    for r in range(1, FUNC_RUNS[paradigm] + 1):
        
        html.write_section_head("Run %d"%r)

        srcdir = os.path.join(ANALYSIS_DIR,paradigm,subj,"model","smoothed","run_%d"%r)
        
        glm_img = os.path.join(srcdir, "design_image.png")
        html.write_text("Design Matrix")
        html.write_image(glm_img)
        html.newline()

        cov_img = os.path.join(srcdir, "stimulus_correlation.png")
        html.write_text("Design Covariance")
        html.write_image(cov_img)
        html.newline()

        res_img = os.path.join(srcdir, "sigmasquareds.png")
        html.write_text("Residual Variance")
        html.write_image(res_img)
        html.newline()

def write_stats_report(html, subj, paradigm):

    for r in range(1, FUNC_RUNS[paradigm] +1):
        html.write_section_head("Run %d"%r)
        srcdir = os.path.join(ANALYSIS_DIR,paradigm,subj,"model","smoothed","stats")
        contrasts = [p.split("/")[-1] for p in glob(os.path.join(srcdir,"*"))]
        contrasts.sort()
        srcdir = os.path.join(ANALYSIS_DIR,paradigm,subj,"model","smoothed","stats")
        for contrast in contrasts:
            html.newline()
            html.write_text(contrast)
            html.write_image(os.path.join(srcdir,contrast,"run_%d"%r,"zstat.png"))
            html.newline()

def write_ffx_report(html, subj, paradigm):
    
    srcdir = os.path.join(ANALYSIS_DIR,paradigm,subj,"fixed_fx", "volume")
    contrasts = [p.split("/")[-1] for p in glob(os.path.join(srcdir,"*"))]
    contrasts.sort()
    html.newline()
    html.write_section_head("Volume")
    for contrast in contrasts:
        html.newline()
        html.write_text(contrast)
        html.write_image(os.path.join(srcdir,contrast,"zstat.png"))
        html.newline()
    srcdir = os.path.join(ANALYSIS_DIR,paradigm,subj,"fixed_fx", "surface")
    contrasts = [p.split("/")[-1] for p in glob(os.path.join(srcdir,"*"))]
    html.newline()
    html.write_section_head("Surface")
    for contrast in contrasts:
        html.newline()
        html.write_text(contrast)
        for hemi in ["lh", "rh"]:
            html.write_images_across([os.path.join(srcdir,contrast,hemi,"%s-%s.png"%(hemi, view)) \
                                        for view in ["lat", "med", "post"]], (300,300))
        html.newline()

def write_homepage():
    
    html = HTMLReport("/mindhive/gablab/fluid/Analysis/Report/fluid_report.html")

    html.write_title("GFluid Imaging Report")
    html.newline()
    html.write_section_head("Subjects")
    ante_subjects = [p.split("/")[-3] for p in glob(os.path.join(DATA_DIR,"gf??","report","report.html"))]
    post_subjects = [p.split("/")[-3] for p in glob(os.path.join(DATA_DIR,"gf??p","report","report.html"))]
    ante_subjects.sort()
    post_subjects.sort()
    for subjlist in ante_subjects, post_subjects:
        html.write_link_table(os.path.join(DATA_DIR,"%s","report","report.html"),"%s",subjlist,10)
        html.newline()

    html.write_section_head("Preprocessing Diagnostics")
    html.newline()
    pars = [p.split("/")[-2] for p in glob(os.path.join(
        REPORT_DIR,"preproc_diagnostics","*","diagnostics.html"))]
    html.write_link_table(os.path.join(REPORT_DIR,"preproc_diagnostics","%s","diagnostics.html"),
                          "%s", pars, 10)


if __name__ == "__main__":
    main()
