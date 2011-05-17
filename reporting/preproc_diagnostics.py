#! /usr/bin/env python
import os
import sys
from os.path import join as pjoin
import shutil
import argparse

import numpy as np
from scipy import stats
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import htmltools


def main(arglist):

    # Set the exclusions thresholds
    # Probably will move to command line?
    motion_thresh = 3.5
    outlier_thresh = 20
    reg_thresh = 0.6

    # Get reference to possible source trees
    data_dir = os.environ["DATA"]
    analysis_dir = os.environ["ANALYSIS"]

    # Get the command line info
    args = parse_args(arglist)

    # Get an array of subjects ids
    if args.subjects and os.path.isfile(args.subjects[0]):
        subjects = np.loadtxt(args.subjects[0], str)
    else:
        subjects = np.array(args.subjects,str)

    if args.run is None:
        run = 1
        paradigm = args.paradigm
    else:
        run = args.run
        paradigm = "%s-run%d"%(args.paradigm, args.run)
        

    # Figure out the focused source and output dirs
    par_dir = pjoin(analysis_dir, "Nipype", args.paradigm)
    report_dir = pjoin(analysis_dir, "Report", "preproc_diagnostics", paradigm)
    try:
        shutil.rmtree(report_dir)
    except:
        pass
    os.makedirs(report_dir)

    # Write a list of the subject considered in these diagnostics
    np.savetxt(pjoin(report_dir, "all_subjects.txt"), subjects, fmt="%s")

    # Plot a histogram of max relative displacement
    write_motion_plot(args.paradigm, par_dir, subjects, run, report_dir)

    # Plot a histogram of number of outlier TRs
    write_outlier_plot(args.paradigm, par_dir, subjects, run, report_dir)

    # Plot a histogram of bbregister costs
    write_reg_plot(args.paradigm, par_dir, subjects, run, report_dir)

    # Get and write list of subjects excluded for max motion 
    maxima = np.array([get_max_motion(s, par_dir, run) for s in subjects])
    motion_subs = subjects[maxima > motion_thresh]
    np.savetxt(pjoin(report_dir, "exclude_motion.txt"), motion_subs, fmt="%s")

    # Get and write a list of subjects excluded for exessive outliers
    outliers = np.array([get_n_outliers(s, par_dir, run) for s in subjects])
    outlier_subs = subjects[outliers > outlier_thresh]
    np.savetxt(pjoin(report_dir, "exclude_outliers.txt"), outlier_subs, fmt="%s")

    # Get and write a list of all subjects to be excluded (the union of motion and outlier exclusions)
    exclude_all = np.array(sorted(list(set(np.concatenate((motion_subs, outlier_subs))))))
    np.savetxt(pjoin(report_dir, "exclude_all.txt"), exclude_all, fmt="%s")

    # Get and write a list of potential problem registrations
    costs = np.array([get_reg_cost(s, par_dir, run) for s in subjects])
    trouble_subs = subjects[costs > reg_thresh]
    np.savetxt(pjoin(report_dir, "high_reg_costs.txt"), trouble_subs, fmt="%s")

    # Write an html summary of the above processing
    write_html_report(args.paradigm, par_dir, subjects, run, motion_thresh, outlier_thresh, report_dir)

def get_max_motion(s, p_dir, run):
    
    return float(open(pjoin(p_dir,s,"preproc/run_%d/max_motion.txt"%run)).readlines()[3].strip())

def get_n_outliers(s, p_dir, run):

    return len(open(pjoin(p_dir,s,"preproc/run_%d/outlier_volumes.txt"%run)).readlines())

def get_reg_cost(s, p_dir, run):
    
    return float(open(pjoin(p_dir,s,"preproc/run_%d/func2anat.mincost"%run)).readline().split()[0])
    
def write_motion_plot(par, p_dir, subjects, run, report_dir):

    maxima = [get_max_motion(s, p_dir, run) for s in subjects]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("%s maximum relative displacement"%par)
    ax.set_xlabel("maximum relative displacement (mm)")
    ax.set_ylabel("number of subjects")
    ax.hist(maxima, 30)

    plt.savefig(pjoin(report_dir, "max_motion_hist.png"))

def write_outlier_plot(par, p_dir, subjects, run, report_dir):

    outliers = [get_n_outliers(s, p_dir, run) for s in subjects]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("%s ART outlier volumes"%par)
    ax.set_xlabel("number of outlier volumes")
    ax.set_ylabel("number of subjects")
    ax.hist(outliers, 30)

    plt.savefig(pjoin(report_dir, "art_outlier_hist.png"))

def write_reg_plot(par, p_dir, subjects, run, report_dir):

    costs = [get_reg_cost(s, p_dir, run) for s in subjects]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("%s bbregister final cost value"%par)
    ax.set_xlabel("final cost value")
    ax.set_ylabel("number of subjects")
    ax.hist(costs, 30)

    plt.savefig(pjoin(report_dir, "bbreg_cost_hist.png"))

def write_html_report(par, p_dir, subjects, run, motion_thresh, outlier_thresh, report_dir):

    html = htmltools.HTMLReport(pjoin(report_dir, "diagnostics.html"))

    html.write_title("%s preprocessing diagnostics"%par)
    html.newline(2)
    html.write_link("subjects involved in preprocessing", "all_subjects.txt")
    html.write_section_head("Motion")
    html.write_image("max_motion_hist.png")
    html.newline()
    maxima = [get_max_motion(s, p_dir, run) for s in subjects]
    maxinfo = stats.describe(maxima)
    html.write_text("min peak motion: %.2f mm"%maxinfo[1][0])
    html.newline()
    html.write_text("max peak motion: %.2f mm"%maxinfo[1][1])
    html.newline()
    html.write_text("mean peak motion: %.2f mm"%maxinfo[2])
    html.newline()
    html.write_text("variance of peak motion: %.2f mm"%maxinfo[3])
    html.newline(2)
    html.write_text("max relative displacement exclusion threshold: %.2f mm"%motion_thresh)
    html.newline()
    try:
        n_exclusions = np.loadtxt(pjoin(report_dir, "exclude_motion.txt"), str).size
    except IOError:
        n_exclusions = 0
    html.write_text("number of subjects excluded for motion: %d"%n_exclusions)
    html.write_link("subjects excluded for motion", "exclude_motion.txt")
    html.newline()
    html.write_section_head("ART Outliers")
    html.newline()
    html.write_image("art_outlier_hist.png")
    html.newline()
    outliers = [get_n_outliers(s, p_dir, run) for s in subjects]
    outinfo = stats.describe(outliers)
    html.write_text("min outliers: %d vols"%outinfo[1][0])
    html.newline()
    html.write_text("max outliers %d vols"%outinfo[1][1])
    html.newline()
    html.write_text("mean outliers %.2f vols"%outinfo[2])
    html.newline()
    html.write_text("variance of outliers: %.2f vols"%outinfo[3])
    html.newline(2)
    html.write_text("art outlier exclusion threshold: %d volumes"%outlier_thresh)
    html.newline()
    try:
        n_exclusions = np.loadtxt(pjoin(report_dir, "exclude_outliers.txt"), str).size
    except IOError:
        n_exclusions = 0
    html.write_text("number of subjects excluded for outliers: %d"%n_exclusions)
    html.write_link("subjects excluded for outliers", "exclude_outliers.txt")
    try:
        n_exclusions = len(np.loadtxt(pjoin(report_dir, "exclude_all.txt"), str))
    except (IOError, TypeError):
        n_exclusions = 0
    html.write_text("number of excluded subjects: %d"%n_exclusions)
    html.newline()
    html.write_link("all automatically excluded subjects", "exclude_all.txt")
    html.newline()
    html.write_section_head("Registration")
    html.write_image("bbreg_cost_hist.png")
    html.newline()
    costs = [get_reg_cost(s, p_dir, run) for s in subjects]
    costinfo = stats.describe(costs)
    html.write_text("min registration cost: %.2f"%costinfo[1][0])
    html.newline()
    html.write_text("max registration cost %.2f"%costinfo[1][1])
    html.newline()
    html.write_text("mean registration cost %.2f"%costinfo[2])
    html.newline()
    html.write_text("variance of registration costs: %.2f"%costinfo[3])
    html.newline(2)
    try:
        n_reg_issues = np.loadtxt(pjoin(report_dir, "high_reg_costs.txt"), str).size
    except IOError:
        n_reg_issues = 0
    html.write_text("number of subjects with potential reg issues: %d"%n_reg_issues)
    html.newline()
    html.write_link("possibly problematic registrations", "high_reg_costs.txt")
    html.newline()
    html.write_link("back to report homepage", "../../fluid_report.html")


def parse_args(arglist):

    parser = argparse.ArgumentParser()
    parser.add_argument("-subjects", nargs="*")
    parser.add_argument("-paradigm")
    parser.add_argument("-run", type=int)

    return parser.parse_args(arglist)

if __name__ == "__main__":
    main(sys.argv[1:])
