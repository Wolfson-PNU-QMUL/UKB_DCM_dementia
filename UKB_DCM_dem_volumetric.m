%%% Sam Ereira (s.ereira@qmul.ac.uk) - last updated 2/3/24

function  UKB_DCM_dem_volumetric(demos)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% FIRST: FIT THE CASE-CONTROL CLASSIFIER TO VOLUMETRIC DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract subsegmental volumes from 40 hippocampal and parahippocampal
% regions
X = UKB_DCM_dem_subseg_volumes(demos);

% Normalise by total intracranial volume
TIV = demos.R_TIV_f25009f25003_2;
X = X./repmat(TIV, 1, size(X,2));

%impute missing data
for i = 1:size(X, 2)
    column_mean = nanmean(X(:, i));
    X(isnan(X(:, i)), i) = column_mean;
end

% Binary outcome (case vs control)
Y = ~strcmp(demos.R_ML_DiagbySess_C42C240Xf41270f20002_Dementia_2,'No');

%Fit K-fold cross-validated predictive model (elastic-net logistic regression)
[metrics] = UKB_DCM_dem_general_classifier(X,Y,10,5);
dt = string(datetime, 'dd-MM-yy_hh:mm:ss');
save(['volumetric_classifier_' dt{1} '.mat'], 'metrics')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SECOND: FIT THE PROGNOSTICATOR TO VOLUMETRIC DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Exclude this subject as unsure when was diagnosed
demos(find(demos.EID == 1595913),:) = [];

%cases only
inds = ~strcmp(demos.R_ML_DiagbySess_C42C240Xf41270f20002_Dementia_2,'No');
demos = demos(inds,:);

% extract subsegmental volumes from 40 hippocampal and parahippocampal
% regions
X = UKB_DCM_dem_subseg_volumes(demos);

% Normalise by total intracranial volume
TIV = demos.R_TIV_f25009f25003_2;
X = X./repmat(TIV, 1, size(X,2));

%impute missing data
for i = 1:size(X, 2)
    column_mean = nanmean(X(:, i));
    X(isnan(X(:, i)), i) = column_mean;
end

% Continuous outcome (time until diagnosis)
Y = str2double(demos.ML_C42C240Xf41270f20002_DementiaAge) - cell2mat(table2cell(demos(:,'R_Age_f21003_2')));

%Fit K-fold cross-validated predictive model (elastic-net linear regression)
[metrics] = UKB_DCM_dem_general_prognosticator(X,Y,10,5);
dt = string(datetime, 'dd-MM-yy_hh:mm:ss');
save(['volumetric_prognosticator_' dt{1} '.mat'], 'metrics')


end