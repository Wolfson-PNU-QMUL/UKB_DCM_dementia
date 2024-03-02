%%% MASTER ANALYSIS PIPELINE for DCM on DMN in dementia cases and controls
%%% from the UK Biobank
%%% Sam Ereira (s.ereira@qmul.ac.uk) - last updated 2/3/24
%%% Code is specific for a specific DMN set of ROIs and the sub-functions
%%% regularly call specific file directories unique to the UK Biobank
%%% dataset for a specific project.  

%
% The preprocessing functions take a maximum framewise displacement (maxFD)
% as an input argument. This is set to 2.4 mm (the voxel size of the UKB
% dataset). Subjects with maximum framewise displacement exceeding this
% value are excluded from the the analysis pipeline. If this number is set
% to a higher value then they won't be excluded. 
%
% DEPENDENCIES
% 1) SPM  - https://www.fil.ion.ucl.ac.uk/spm/software/download/
% 2) Glmnet - https://hastie.su.domains/glmnet_matlab/

clear
clc

%user must specify path and directory for functional MRI data
fMRI_datadir = '';

%user must specify path and directory for structural MRI data
sMRI_datadir = '';

%user must specify path for SPM12 (software)
SPM_dir = '';

%Load subject IDs and demographic data
incident = readtable('UKBW_rsfmri_F_DemIncident_F10_280623.csv');
prevalent = readtable('UKBW_rsfmri_F_DemPrevalent_F10_280623.csv');

demos = [incident;prevalent];
subjects = [incident.EID;prevalent.EID];

%Maximum acceptable framewise displacement is 2.4 mm (voxel size of this dataset)
maxFD = 2.4;

%Specify which analysis step you want to run (only does one at at time)
job = input(sprintf("Select a pipeline step to implement:\n1 = Preprocess Structurals\n2 = Preprocess Functionals\n3 = Extract Timeseries " + ...
    "\n4 = Fit 1st level DCMs\n5 = Fit 2nd level DCM (PEB) as function of dementia incidence and classify \n6 = Fit 2nd level DCM (PEB) as a function of time until diagnosis (cases only) and predict" + ...
    "\n7 = Analysis on modifiable risk factors\n8 = Volumetric comparator analysis\n9 = Functional connectivity comparator analysis\n"));

switch job
    case 1
        UKB_DCM_dem_structural_MRI_preprocess(subjects, sMRI_datadir, SPM_dir) 
    case 2
        UKB_DCM_dem_functional_MRI_preprocess(subjects, fMRI_datadir, sMRI_datadir, maxFD) 
    case 3
        UKB_DCM_dem_extract_timeseries(subjects, fMRI_datadir, maxFD) 
    case 4
        UKB_DCM_dem_firstlevelDCM(subjects, fMRI_datadir, maxFD) 
    case 5
        UKB_DCM_dem_EC_classifier(subjects, fMRI_datadir, demos) 
    case 6
        UKB_DCM_dem_EC_prognosticator(subjects, fMRI_datadir, demos) 
    case 7
        preclin_only = 0; %switch to 1 to exclude those who had dementia at time of scan
        UKB_DCM_dem_ModRiskFactors(preclin_only, incident, prevalent, SPM_dir) 
    case 8
        UKB_DCM_dem_volumetric(demos)
    case 9
        UKB_DCM_dem_functional_connectivity(demos)
end