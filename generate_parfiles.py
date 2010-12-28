#! /usr/bin/env python
import os
import sys
from glob import glob
import numpy as np
import scipy.io as scio

def main():

    # Get subject from first command line arg
    try:
        subject = sys.argv[1]
    # Print usage and quit if we didn't get a subject
    except IndexError:
        sys.exit("\tUSAGE: generate_parfiles.py subject_id")

    # Figure out ante/post based on subject id
    if subject.endswith("p"):
        day = 2
    else:
        day = 1
    

    matfiletemp = dict(nback="sub_%s_day%d_run%d_DUALNBACK.mat",
                       rt="sub_%s_r_%d_TEST_RT_FMRI.mat",
                       iq="sub_%s_r_?_s_scan%d_BLOCK_CATTELL.mat",
                       mot="%s_Response*.mat")
    
    # Get directory locations
    datadir = "/mindhive/gablab/fluid/Data/%s/"%subject
    if not os.path.isdir(datadir):
        sys.exit("\nError: Data directory for subject %s does not exist\n"%subject)
    srcdir = os.path.join(datadir, "fmri_bhvl")
    trgdir = os.path.join(datadir, "parfiles")
    # Make the parfile directory if it doesn't exist
    if not os.path.exists(trgdir): 
        os.mkdir(trgdir)

    gen_nback(matfiletemp["nback"], day, subject, srcdir, trgdir)

    gen_iq(matfiletemp["iq"], day, subject, srcdir, trgdir)
    
    gen_mot(matfiletemp["mot"], day, subject, srcdir, trgdir)

    #gen_rt(matfiletemp["rt"], subject, srcdir, trgdir)

def gen_nback(matfiletemplate, day, subject, srcdir, trgdir):
    print "\nNBack\n======="
    for run in range(1, 5):
        try:
            matfile = os.path.join(srcdir, matfiletemplate%(subject, day, run))
            print "Reading %s"%matfile
            stE = scio.loadmat(matfile,struct_as_record=False,squeeze_me=True)["stE"]

            nblocks = len(stE.Block)
            trialtypes = {0:"zero",1:"easy",2:"medium",3:"hard"}
            for nback in [0, 1, 2, 3]:
                parfile = os.path.join(
                    trgdir, "NBack_r%d_d%d_%s_%s.txt"%(
                        run,day,trialtypes[nback],subject))
                print "   Writing ", parfile
                fid = open(parfile, "w")
                for block in range(nblocks):
                    if stE.Block[block].iNumBack == nback:
                        onset = stE.Block[block].Trials.dActStimOnset[0]
                        duration = np.subtract(stE.Block[block].dActFeedbackStart, onset)
                        fid.write("%.2f\t%.2f\t%d\n"%(onset,duration,1))
                fid.close()
            instructfile = os.path.join(
                trgdir, "NBack_r%d_d%d_inst_%s.txt"%(run,day,subject))
            fid = open(instructfile,"w")
            print "   Writing %s"%instructfile
            for block in range(nblocks):
                fid.write("%.2f\t%.2f\t%d\n"%(stE.Block[block].dActInstStart,3,1))
            fid.close()

        # Catch the error from a missing .mat file
        except IOError:
            if not os.path.exists(matfile):
                print "ERROR: could not read %s"%matfile
            else:
                raise

def gen_mot(matfiletemplate, day, subject, srcdir, trgdir):
    print "\nMOT\n======="
    for matfile in glob(os.path.join(srcdir, matfiletemplate%subject)):
        print "Reading %s"%matfile
        try:
            D = scio.loadmat(matfile, struct_as_record=False, squeeze_me=True)["D"]
            runtype = {129:"Block",193:"Jitter"}[len(D.actualEventStartTime)]
            for speed in range(1,5):
                parfile = os.path.join(trgdir, "MOT_%s_r1_d1_speed%d_%s.txt"%(runtype,speed,subject))
                print "   Writing %s"%parfile
                fid = open(parfile,"w")
                for i, time in enumerate(D.actualEventStartTime):
                    if D.eventCond[i] == speed:
                        if D.eventType[i] == 3:
                            now = D.actualEventStartTime[i]
                            then = D.actualEventStartTime[i+2]
                            fid.write("%.2f\t%.2f\t1\n"%(now, then-now))
                fid.close()
            parfile = os.path.join(trgdir, "MOT_%s_r1_d1_resp_%s.txt"%(runtype,subject))
            print "   Writing %s"%parfile
            fid = open(parfile,"w")
            for i, time in enumerate(D.actualEventStartTime):
                if D.eventType[i] == 5:
                    now = D.actualEventStartTime[i]
                    then = D.actualEventStartTime[i+1]
                    fid.write("%.2f\t%.2f\t1\n"%(now, then-now))
            fid.close()

        # Catch the error from the missing .mat file
        except IOError:
            if not os.path.exists(matfile):
                print "ERROR: could not read %s"%matfile
            else:
                raise

def gen_iq(matfiletemplate, day, subject, srcdir, trgdir):
    print "\nIQ\n======="
    try:
        matfile = glob(os.path.join(srcdir, matfiletemplate%(subject, day)))[0]
        print "Reading %s"%matfile
        stE = scio.loadmat(matfile, struct_as_record=False, squeeze_me=True)["stE"]
        diffdict = dict(easy="LOW",hard="HIGH")
        for diff in diffdict:
            parfile = os.path.join(
                trgdir, "IQ_r%d_d%d_%s_%s.txt"%(1,day,diff,subject))
            print "   Writing ", parfile
            fid = open(parfile, "w")
            for i, block in enumerate(stE.Block.sBlockType):
                if block == diffdict[diff]:
                    fid.write("%.3f\t%d\t%d\n"%(stE.Block.dActBlockOnset[i],45,1))
            fid.close()

    # Catch the error from a missing .mat file
    except IOError:
        if not os.path.exists(matfile):
            print "ERROR: could not read %s"%matfile
        else:
            raise

def gen_rt(matfiletemplate, subject, srcdir, trgdir):
    print "\nRT\n======="
    for run in range(1,3):
        try:
            # Read the mat file
            matfile = os.path.join(srcdir, matfiletemplate%(subject, run))
            print "Reading %s"%matfile
            stE = scio.loadmat(matfile, struct_as_record=False, squeeze_me=True)["stE"]

            # If the length of time between trials is > 4, assume a new block
            deltas = np.subtract(stE.Trials.dActInitFix[1:], stE.Trials.dActInitFix[:-1])
            # Make the first delta very large as we know it's a new block
            deltas = np.hstack((np.array(999),deltas))
            # Find the indices for the first fixations of new blocks
            id = np.argwhere(deltas > 4)
            # Get the vector of onsets
            onsets = stE.Trials.dActInitFix[id]
            # Get the vector of trials types
            blocktype = stE.Trials.iTrialType[id]

            # Make the parfiles for events of interest
            typedict = {1:"simple",2:"choice"}
            for trialtype in [1, 2]:
                parfile = os.path.join(
                    trgdir, "RT_r%d_%s_%s.txt"%(run,typedict[trialtype],subject))
                print "   Writing ", parfile
                fid = open(parfile, "w")
                for i, onset in enumerate(onsets):
                    if blocktype[i] == trialtype:
                        fid.write("%.2f\t%.2f\t%d\n"%(onset, 30, 1))
                fid.close()

            # Make instructions parfile
            parfile = os.path.join(trgdir, "RT_r%d_%s_%s.txt"%(run,"inst", subject))
            print "   Writing ", parfile
            fid = open(parfile, "w")
            for onset in stE.dActBlockInstStart:
                fid.write("%.2f\t%.2f\t%d\n"%(onset,2,1))
            fid.close()

        except IOError:
            if not os.path.exists(matfile):
                print "ERROR: could not read %s"%matfile
            else:
                raise error

if __name__ == "__main__":

    main()
