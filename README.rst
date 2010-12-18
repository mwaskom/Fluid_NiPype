Fluid Intelligence Project Analysis Code
========================================

Functional Analysis
-------------------

``fluid_fmri.py`` -- 
Main interface for functional analysis. 
Uses code from workflows and experiment packages.

``fluid_resting.py`` -- 
Interface for resting state preprocessing and registration. 
Uses code from workflows package.


Workflows
---------
Package with workflow modules.
*Should* be entirely independent of anything Gfluid relevant.

``preproc.py`` -- 
Preprocessing.

``fsl_model.py`` -- 
First-level FSL modelfitting .

``spm_model.py`` -- 
First-level SPM modelfitting.

``fixed_fx.py`` -- 
Second-level (within-subject) fixed-effect analysis.

``resting_preproc.py`` -- 
Preprocessing for resting-state scans.

``vbm_preproc.py``
Preprocssing for arbitrary VBM style analysis.

Experiment
----------
Package with modules containing paradigm information.

``iq.py`` -- 
fMRI paradigm with psychometric tasks.

``nback.py`` -- 
Dual (auditory letter-matching/spatial) N-Back task.

``mot_block.py`` -- 
Blocked multiple object tracking.

``mot_jitter.py`` -- 
Event-related multiple object tracking.

Structural Analysis
-------------------

``fluid_flash.py`` -- 
Analysis of multispectral FLASH images.

``fluid_dti.py`` -- 
Analysis of DTI-based images.

Misc
----

``unpack_dicoms.py`` -- 
NiPype script to convert DICOM files to nifti/mgh format and perform 
structural preprocessing.

``generate_parfiles.py`` -- 
Script to read Psychtoolbox experiment files and create FSL parfiles.

``fluid_utility.py`` -- 
Assorted functions used in our NiPype scripts.

Note
----

The code contained in this repository is subject to change at any time.
(Although I would like to think that the processing and registration workflows 
are pretty stable).
