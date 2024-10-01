% Code written by Sam Ereira (s.ereira@qmul.ac.uk) - last updated 2/3/24
% Adapted from Thomas et al. 2023  - github.com/gecthomas/Spectral_DCM_in_PD_VH
%
% Fits a second-level (i.e. group level) DCM model and models differences amongst cases
% as a function of time until diagnosis and then does bayesian model
% reduction/averaging. The result of this second-level model is a connectivity parameter for
% each DCM connection for each between-subject variable. The DCM here has a
% total of 100 connectivity parameters in the A matrix (fully connected
% network of 10 DMN nodes). Therefore the final number of parameters that
% we get from this second-level model will be 100 multipled by however many
% between-subjects variables we include in the between-subjects design
% matrix. Here we include 5 columns in the between-subjects design matrix
% and there are therefore 500 parameters estimated.
%
% Then applies elastic-net regularised linear
% regression using K-fold cross validation to see if time until diagnosis
% can be predicted from effective connectivity
%

function UKB_DCM_dem_EC_prognosticator(subjects, fMRI_datadir,  demos)

K = 10; %10 folds of outer cross validation in K-fold
k = 5; %5 folds of inner cross validation

%Mean, Age at diagnosis, age of scan, sex
%We want to transform age at diagnosis into time until diagnosis
labels = {'Ctrls','ML_C42C240Xf41270f20002_DementiaAge', 'R_Age_f21003_2', 'R_Sex_f31'};

%covariates as they'll be labelled for visualisation
labels_vis = {'Mean','TimeUntilDiagnosis', 'Age', 'Sex', 'FramewiseDisplacement'};

%same for all subjects
funct_ID = '_20227_2_0';


%%%%%%% construct GCM from individual DCMs %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Quality check first. Only include DCMs that meet fit criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
excluded = [];
GCM = []; GCM_sj_list = [];
j = 0; %counter for converged DCMs
for i = 1:length(subjects)
    converged = [];
    EID = subjects(i);
    if EID == 1595913 %UNSURE OF DEMENTIA AGE SO EXCLUDE THIS SUBJECT
        continue
    end
    funct_data_path = [fMRI_datadir num2str(EID) funct_ID];

    %Don't include any controls. Only cases
    if ~strcmp(demos.R_ML_DiagbySess_C42C240Xf41270f20002_Dementia_2(find(sum(demos.EID == EID,2))), 'No')

        %screen DCM before putting to 2nd level
        load([funct_data_path '/fMRI/Full_DCM/Diagnostics'], 'converged');

        if converged == 1
            j = j + 1;
            load([funct_data_path '/fMRI/Full_DCM/DCM_full_estim']);
            GCM{j,1} = DCM{1};
            GCM_sj_list(j,1) = EID;
        else
            excluded = [excluded, i];
        end
    end
end

fprintf(['Fitting PEB to ' num2str(j) ' out of a possible total of ' num2str(i) ' DCMs'])

%create between-subjects design matrix to model group level effects
N = length(GCM_sj_list); %number of subjects

%number of between-subjects variable
M = length(labels);

%initialise design
X = nan(N, M);

%First column is ones (to represent baseline, i.e. controls)
X(:, 1) = ones(N,1);

%indices of chosen subjects within demographics table for the full cohort
[~,sj_inds] = intersect(demos.EID, GCM_sj_list, 'stable');

for col = 2:M %go through each covariate in turn
    data = (table2cell(demos(sj_inds,labels{col})));
    if strcmp(labels{col}, 'ML_C42C240Xf41270f20002_DementiaAge')
        %transform this from age at diagnosis to time until diagnosis
        config_data =  str2double(data) - cell2mat(table2cell(demos(sj_inds,'R_Age_f21003_2')));

        %mean-centre
        X(:, col) = config_data - mean(config_data);
    else
        %mean-centre
        X(:, col) = cell2mat(data) - mean(cell2mat(data)); %mean centre
    end
end

% Add 5th and final column to X (framewise displacement) to model motion
tmp = [];
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
[PEB, RCM] = spm_dcm_peb(GCM,M,{'A'}); %RCM is a 'reduced'/'refined' array of individual DCMs which have been constrained by group-level priors

%Bayesian model reduction and average using automatic greedy search over parameters
%Prunes connections
BMA = spm_dcm_peb_bmc(PEB);

dt = string(datetime, 'dd-MM-yy_hh:mm:ss');
file = ['Prognosticator_DCM_BMA_' dt{1} '.mat'];

fprintf('PEB/BMA complete. Saving files...')
save(file,'BMA', 'PEB', 'RCM', 'GCM','-v7.3')


%%%%%%%%%%%%% Now see if we can actually use the surviving connections from
%%%%%%%%%%%%% the Bayesian model reduction to predict date of dementia
%%%%%%%%%%%%% diagnosis using leave-one-out cross-validation
clearvars -except BMA PEB GCM RCM demos sj_inds K k
clc
close all
fprintf('Attempting to predict date of diagnosis using effective connectivity...')

%number of regions
nR = length(RCM{1,1}.Y.name);

%Only use connections with very strong posterior evidence of being non-zero
inds = find(BMA.Pp(((nR^2)+1) :(2*(nR^2))  )>0.99);
[to,from] = ind2sub([nR nR],inds);

Y = []; X = [];
for isj = 1:length(GCM)
    for j = 1:length(to)

        %Data features (DCM parameters)
        X(isj, j) = RCM{isj}.Ep.A(to(j), from(j));
        
        %N.B. For future users
        %This analysis makes use of RCM (parameters refined by group-level priors) which may bias subsequent results
        %One alternative would have bene to exclude the primary covariate of interest (time until diagnosis) from the original PEB analysis
        %However, the most appropriate analysis here is probably to use Bayesian LOO-CV within the PEB framework, rather than using 
        %classical inference, outside of the PEB framework, as we have done here

        %Response variable to predict (time until diagnosis)
        Y(isj,1) =  str2double(cell2mat(table2cell(demos(sj_inds(isj), 'ML_C42C240Xf41270f20002_DementiaAge')))) - cell2mat(table2cell(demos(sj_inds(isj),'R_Age_f21003_2')));

    end
end

%Fit cross-validated predictive model (elastic-net linear regression)
[metrics] = UKB_DCM_dem_general_prognosticator(X,Y,K,k);

dt = string(datetime, 'dd-MM-yy_hh:mm:ss');
save(['EC_prognosticator_' dt{1} '.mat'], 'metrics')

