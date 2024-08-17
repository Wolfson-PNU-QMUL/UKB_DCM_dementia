% Code written by Sam Ereira (s.ereira@qmul.ac.uk) - last updated 2/3/24
% Adapted from Thomas et al. 2023  - github.com/gecthomas/Spectral_DCM_in_PD_VH
%
% Fits a second-level (i.e. group level) DCM model and models differences
% between cases and controls and then does bayesian model
% reduction/averaging. The result of this second-level model is a connectivity parameter for
% each DCM connection for each between-subject variable. The DCM here has a
% total of 100 connectivity parameters in the A matrix (fully connected
% network of 10 DMN nodes). Therefore the final number of parameters that
% we get from this second-level model will be 100 multipled by however many
% between-subjects variables we include in the between-subjects design
% matrix. Here we include 5 columns in the between-subjects design matrix
% and there are therefore 500 parameters estimated.
% 
% Then applies elastic-net regularised logistic
% regression using K-fold cross validation to see if dementia
% incidence can be predicted from effective connectivity
%
%
% NOTE ON THE BETWEEN SUBJECTS DESING MATRIX AND THE IMPORTANCE OF MEAN
% CENTERING
% The between-subjects design matrix here does NOT mean centre the primary
% variable of interest (case vs. control).
% That means that the first column in the design matrix (a column of 1s)
% represents the average connectivity matrix in control participants and
% the second column in the design matrix represents the DIFFERENCE in
% connectivity between cases and controls. The average connectivity of
% cases can be derived by summing the outputted connectivity matrices for the first
% and second columns of the between-subjects design matrix
%
% If we HAD decided to mean centre the second column of the design matrix
% then the first column would represent the average connectivity matrix
% (commonalities) across ALL subjects (cases and controls). The second column would then
% represent how cases and controls both deviate from this overall group
% average. Average connectivity for controls would be derived by
% subtracting the outputted connectivity matrix for the 2nd design-matrix
% column from that of the 1st design-matrix column. The average
% connectivity for cases would be derived by ADDING these two connectivity
% matrices.

function UKB_DCM_dem_EC_classifier(subjects, fMRI_datadir, demos)

K = 10; %10 folds of outer cross validation in K-fold
k = 5; %5 folds of inner cross validation

%between-subjects variables in design matrix: dementia status, age, sex +
%Framewise displacement
labels = {'Ctrls', 'R_ML_DiagbySess_C42C240Xf41270f20002_Dementia_2','R_Age_f21003_2', 'R_Sex_f31'};

%covariates as they'll be labelled for visualisation
labels_vis = {'Controls', 'Cases - Controls', 'Age', 'Sex', 'Framewise Displacement'};

funct_ID = '_20227_2_0'; %same for all subjects


%%%%%%% construct GCM from individual DCMs %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Quality check first. Only include DCMs that satisfactorily converged
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
excluded = [];
GCM = [];
j = 0; %counter for converged DCMs

for i = 1:length(subjects)
    converged = [];
    EID = subjects(i);
    funct_data_path = ['/data/Wolfson-PNU-dementia/UKB/Imaging_rsfMRI_UKB_DCM/rfMRI_CCDem_2/' num2str(EID) funct_ID];


    %screen DCM before putting to 2nd level
    load([funct_data_path '/fMRI/Full_DCM/Diagnostics'], 'converged');

    if converged == 1
        j = j + 1;
        load([funct_data_path '/fMRI/Full_DCM/DCM_full_estim']);
        GCM{j,1} = DCM{1};
        GCM_sj_list(j,1) = EID; %ID

        %case (1) or control (0) or prevalent (2)
        GCM_sj_list(j,2) = ~strcmp(demos.R_ML_DiagbySess_C42C240Xf41270f20002_Dementia_2(i), 'No');
        if strcmp(demos.R_ML_DiagbySess_C42C240Xf41270f20002_Dementia_2(i), 'Yes')==1
            GCM_sj_list(j,2) = 2; %label for prevalent cases
        end

        %position in array for indexing later
        GCM_sj_list(j,3) = i;

    else
        excluded = [excluded, i];
    end
end

fprintf(['Fitting PEB to ' num2str(j) ' out of a possible total of ' num2str(i) ' DCMs'])

%create between-subjects design matrix to model group level effects
N = length(GCM_sj_list); %number of subjects
M = length(labels); %number of between-subjects variables
X = nan(N, M); %initialise design matrix

%First column is ones (to represent baseline, i.e. controls)
X(:, 1) = ones(N,1); 

%indices of chosen subjects within demographics table for the full cohort
[~,sj_inds] = intersect(demos.EID, GCM_sj_list, 'stable');

counter = 0;
for col = 2:M %go through each between-subjects variable in turn

    data = table2cell(demos(sj_inds,labels{col}));

    if strcmp(labels{col}, 'R_ML_DiagbySess_C42C240Xf41270f20002_Dementia_2')
        counter = counter + 1;
        config_data = ones(length(data), 1); %default to case
        config_data(find(strcmp(data,'No'))) =0; %select controls
        X(:, col) = config_data; 
    else
        %mean-centre other columns
        X(:, col) = cell2mat(data) - mean(cell2mat(data)); 
    end
end

% Add 5th and final column to X (framewise displacement) to model motion
for isj = 1:length(sj_inds)
    EID = table2array(demos(sj_inds(isj), 1));
    funct_data_path = [fMRI_datadir num2str(EID) funct_ID];
    regressors_directory = [funct_data_path '/fMRI/regressors'];
    load([regressors_directory '/Framewise_Displacement.mat'], 'FD');
    tmp(isj,1) = mean(FD);
    clear FD
end
X = [X, zscore(tmp)];

% PEB global settings
M = struct();
% M.Q is the choice of precision components to use. By setting this to
% single, there will be shared variance parameter between all connectivity
% parameters. If it is set to 'all' then there is a separate covariance
% component for every connectivity parameter. This makes the model very
% complex with a huge number of parmaeters and takes far too long to
% estimate. For simplicity we choose 'single'
M.Q      = 'single';
M.maxit  = 512;
M.X      = X;
M.Xnames = labels_vis;

%fit PEB (2nd level DCM model with all connections)
PEB = spm_dcm_peb(GCM,M,{'A'});

%Bayesian model reduction and average using automatic greedy search over parameters
%Prunes connections
BMA = spm_dcm_peb_bmc(PEB);

dt = string(datetime, 'dd-MM-yy_hh:mm:ss');
file = ['Classifier_DCM_BMA_' dt{1} '.mat'];

fprintf('PEB/BMA complete. Saving files...')
save(file,'BMA', 'PEB', 'GCM', 'GCM_sj_list','-v7.3')


%%%%%%%%%%%%% Now see if we can actually use the surviving connections from
%%%%%%%%%%%%% the Bayesian model reduction to predict dementia status using
%%%%%%%%%%%%% K-fold cross-validation
clearvars -except BMA PEB GCM GCM_sj_list K k
clc
close all
fprintf('Attempting to predict dementia status using effective connectivity...')

%number of regions
nR = length(GCM{1,1}.Y.name);

%Only use connections with very strong posterior evidence of being non-zero
inds = find(BMA.Pp(((nR^2)+1) :(2*(nR^2))  )>0.99);
[to,from] = ind2sub([nR nR],inds);

Y = []; X = [];
for isj = 1:length(GCM)
    for j = 1:length(to)

        %Data features (DCM parameters)
        X(isj, j) = GCM{isj}.Ep.A(to(j), from(j));

        %Response variable to predict (case vs. control)
        Y(isj,1) = PEB.M.X(isj,2); 
    end
end

%Fit cross-validated predictive model (elastic-net logistic regression)
[metrics] = UKB_DCM_dem_general_classifier(X,Y,K,k, GCM_sj_list);

dt = string(datetime, 'dd-MM-yy_hh:mm:ss');
save(['EC_classifier_' dt{1} '.mat'], 'metrics')

end
