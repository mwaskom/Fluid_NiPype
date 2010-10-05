Fluid Intelligence Project Analysis Code
========================================

Functional Analysis
-------------------

``fluid_nipype.py:`` -- 
Main interface for functional analysis

``fluid_source.py`` -- 
Source utility nodes

``fluid_preproc.py`` -- 
Preprocessing

``fluid_fsl_model.py`` -- 
First-level FSL modelfitting 

``fluid_spm_model.py`` -- 
First-level SPM modelfitting

``fluid_fixed_fx.py`` -- 
Second-level (within-subject) fixed-effect analysis

Experiment Modules
------------------

``iq_experiment.py`` -- 
fMRI paradigm with psychometric tasks

``nback_experiment.py`` -- 
Dual (auditory letter-matching/spatial) N-Back task

``mot_block_experiment.py`` -- 
Blocked multiple object tracking

``mot_jitter_experiment.py`` -- 
Event-related multiple object tracking 

``rt_exeriment.py`` -- 
Simple and choice reaction time tasks

Structural Analysis
-------------------

``fluid_flash.py`` -- 
Analysis of multispectral FLASH images

Misc
----

``unpack_dicoms.py`` -- 
NiPype script to convert DICOM files to nifti/mgh format

``fetch_fluid_dicoms.py`` -- 
Wrapper for ``fetch_dicoms`` script that handles our directory structure

``generate_parfiles.py`` -- 
Script to read Psychtoolbox experiment files and create FSL parfiles

``fluid_utility_funcs.py`` -- 
Assorted functions used in our NiPype scripts

Note
----

The code contained in this repository is subject to change at any time.

