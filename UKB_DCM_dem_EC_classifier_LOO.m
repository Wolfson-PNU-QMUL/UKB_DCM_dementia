% Code written by Sam Ereira (s.ereira@qmul.ac.uk) - last updated 2/3/24
% Adapted from Thomas et al. 2023  - github.com/gecthomas/Spectral_DCM_in_PD_VH
% Same as UKB_DCM_dem_EC_classifier but uses leave-out-out (LOO) instead of
% K-fold. This is used for generating an effective connectivity (EC) index
% for each subject.

function UKB_DCM_dem_EC_classifier_LOO(subjects, fMRI_datadir, demos)


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

%option to exclude prevalent cases and their matched controls
if preclin_only== 1
    subjects = subjects(1:891)
    demos = demos(1:891, :);
end

for i = 1:length(subjects)
    converged = [];
    EID = subjects(i);
    funct_data_path = [fMRI_datadir num2str(EID) funct_ID];

    %screen DCM before putting to 2nd level
    load([funct_data_path '/fMRI/Full_DCM/Diagnostics'], 'converged');

    if converged == 1
        j = j + 1;
        load([funct_data_path '/fMRI/Full_DCM/DCM_full_estim']);
        GCM{j,1} = DCM{1};
        GCM_sj_list(j,1) = EID
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
[PEB,RCM] = spm_dcm_peb(GCM,M,{'A'});

%Bayesian model reduction and average using automatic greedy search over parameters
%Prunes connections
BMA = spm_dcm_peb_bmc(PEB);

dt = string(datetime, 'dd-MM-yy_hh:mm:ss');
file = ['Classifier_DCM_BMA_PreClinOnly_' num2str(preclin_only) '_' dt{1} '.mat'];

fprintf('PEB/BMA complete. Saving files...')
save(file,'BMA', 'PEB', 'RCM', 'GCM_sj_list','-v7.3')


%%%%%%%%%%%%% Now see if we can actually use the surviving connections from
%%%%%%%%%%%%% the Bayesian model reduction to predict dementia status using
%%%%%%%%%%%%% leave-one-out cross-validation
clearvars -except BMA PEB RCM preclin_only
clc
close all
fprintf('Attempting to predict dementia status using effective connectivity...')

%number of regions
nR = length(RCM{1,1}.Y.name);

%Only use connections with very strong posterior evidence of being non-zero
inds = find(BMA.Pp(((nR^2)+1) :(2*(nR^2))  )>0.99);
[to,from] = ind2sub([nR nR],inds);

Y = []; X = [];
for isj = 1:length(RCM)
    for j = 1:length(to)

        %Data features (DCM parameters)
        X(isj, j) = RCM{isj}.Ep.A(to(j), from(j));

        %Response variable to predict (case vs. control)
        Y(isj,1) = PEB.M.X(isj,2);
    end
end

pred = [];opt_lam = []; opt_alpha = []; lambdas = [];
alpha_list = 0:0.2:1; %list of alpha (mixing) parameters for elastic net
pars = struct('fdev',-1); %ensures that we use the full sequence of 100 lambdas each time
glmnetControl(pars);
weights = ones(length(RCM),1);
weights(Y==1) = 1-sum(Y==1)./numel(Y); weights(Y==0) = 1-sum(Y==0)./numel(Y);

for isj = 1:length(RCM)
    tic
    fprintf('Elastic net logistic regression, cross-val fold: %s/%s...  ', num2str(isj), num2str(length(RCM)))

    %reset to default options to use full lambda sequence
    options = glmnetSet;

    %Define test-set as one single left-out subject
    Xtest = X(isj, :);

    %Define training set as remainder of subjects
    Xtrain = X; Xtrain(isj,:) = [];
    Ytrain = Y; Ytrain(isj) = [];
    w_outer = weights; w_outer(isj) = [];

    prob_tmp = nan(size(Ytrain,1), numel(alpha_list), 100);lambdas = prob_tmp;
    for jsj = 1:size(Ytrain,1) %Inner loop of further LOO cross-validation
        for alpha_ind = 1:numel(alpha_list) %Loop over hyperparameters to optimise
            options.alpha = alpha_list(alpha_ind);

            %Define inner test set (validation set) as one left out subject
            %from the outer train set
            X_val_test = Xtrain(jsj, :);

            %Define inner train set as remainder of train set
            X_val_train = Xtrain; X_val_train(jsj,:) = [];
            Y_val_train = Ytrain; Y_val_train(jsj) = [];
            w_inner = w_outer; w_inner(jsj) = [];

            %Train logistic regression model
            options.weights = w_inner;
            fit = glmnet(X_val_train, Y_val_train, 'binomial', options);

            %add padding if glmnet stopped short of 100 lambdas
            offset = 100 - size(fit.beta,2);
            if offset >0
                tmp = fit.beta(:,1); pad = repmat(tmp, 1, offset); fit.beta = [pad, fit.beta];
                tmp = fit.lambda(1); pad = repmat(tmp, offset, 1); fit.lambda = [pad;fit.lambda];
                tmp = fit.a0(1); pad = repmat(tmp, 1, offset); fit.a0 = [pad, fit.a0];
            end

            %Output probabilities from regression models...P(case)
            prob_tmp(jsj, alpha_ind,:) = glmnetPredict(fit, X_val_test,[],'response');

            %save lambda sequences as they are always slightly different
            lambdas(jsj,alpha_ind,:) = fit.lambda;
        end
    end


    % Fit a ROC curve for each combination of hyperparameters
    x_tmp = cell(100, numel(alpha_list)); y_tmp = cell(100, numel(alpha_list));
    T_tmp = cell(100, numel(alpha_list)); AUC_tmp = zeros(100, numel(alpha_list));


    for l = 1:size(prob_tmp, 3)
        for alpha = 1:numel(alpha_list)
            [x_tmp{l,alpha},y_tmp{l,alpha},T_tmp{l,alpha},AUC_tmp(l,alpha)] = perfcurve(Ytrain,prob_tmp(:,alpha,l),1);
        end
    end

    %find optimal hyperparameters by maximising AUC-ROC
    [~,ind] = max(AUC_tmp(:));
    [opt_lam_ind, opt_alpha_ind] = ind2sub([size(AUC_tmp,1), size(AUC_tmp,2)],ind);
    opt_lam(isj) = mean(lambdas(:,opt_alpha_ind,opt_lam_ind));
    opt_alpha(isj) = alpha_list(opt_alpha_ind);

    %save a ROC curve for each outer fold of LOO CV
    x(:, isj) = x_tmp{opt_lam_ind, opt_alpha_ind}; y(:, isj) = y_tmp{opt_lam_ind, opt_alpha_ind};

    %find optimal decision boundary for each outer fold of LOO CV
    [~,ind] = max(y(:,isj)-x(:,isj));
    opt_boundary(isj) = T_tmp{opt_lam_ind, opt_alpha_ind}(ind);

    %train classifier on full training set using these optimised parameters
    options.lambda = opt_lam(isj); options.alpha = opt_alpha(isj);
    options.weights = w_outer;
    fit = glmnet(Xtrain, Ytrain, 'binomial', options);

    %test classifier on left out test set and save p(case) and binary prediction
    prob(isj) = glmnetPredict(fit, Xtest,[],'response');
    pred(isj) = prob(isj)>opt_boundary(isj);

    fprintf(['\n' num2str(toc) ' seconds elapsed\n'])
    fprintf(['\n optimal alpha is ' num2str(opt_alpha(isj)) '\n'])
    fprintf(['\n max AUC is ' num2str(max(AUC_tmp(:))) '\n'])
end

%%%% Generate final model performance metrics %%%%

% AUC, sensitivity, specificity and predictive values
[metrics.AUC.x,metrics.AUC.y,metrics.AUC.T,metrics.AUC.AUC] = perfcurve(Y,prob,1);
metrics.sensitivity = sum(pred(Y==1)==1)/sum(Y==1);
metrics.specificity = sum(pred(Y==0)==0)/sum(Y==0);
metrics.PPV = sum(Y(pred==1)==1)/sum(pred==1);
metrics.NPV = sum(Y(pred==0)==0)/sum(pred==0);

% Point biserial Correlation
[metrics.corr.r,~,metrics.corr.p,metrics.corr.ci] = pointbiserial(Y,prob);

save(['LOOCV_EC_classifier_weighted_PreClinOnly_' num2str(preclin_only) '.mat'], 'metrics','x', 'y', 'opt_boundary','opt_lam', 'opt_alpha', 'prob', 'pred', 'X', 'Y')


end
