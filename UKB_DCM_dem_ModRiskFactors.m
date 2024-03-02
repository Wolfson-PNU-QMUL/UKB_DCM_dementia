function UKB_DCM_dem_ModRiskFactors(preclin_only, incident, prevalent, SPM_dir)

%Remove path for SPM directory as this function otherwise fitglm doesn't
%work properly
rmpath(SPM_dir)

%Load effective connectivity index (i.e. classifier-output of P(dementia)
%from a leave-one-out cross-validated classifier 
% (pre-generated using UKB_DCM_dem_EC_classifier_LOO.m) 
if preclin_only == 1
    load('LOOCV_EC_classifier_weighted_PreClinOnly_1.mat', 'prob')
    EC_index = prob';
    demos = incident; %load demographic data
else
    load('LOOCV_EC_classifier_weighted_PreClinOnly_0.mat', 'prob')
    EC_index = prob';
    demos = [incident;prevalent]; %load demographic data
end


%initialise variables for social-isolation composite variable
X = [];

%visits from friends/family. Convert into number of visits per year
tmp = demos.R_SI_FreqVisits_f1031_2;
config = nan(length(tmp), 1);
config(contains(tmp, 'Never')) = 0;
config(contains(tmp, 'few months')) = 3;
config(contains(tmp, '1x/month')) = 12;
config(contains(tmp, '1x/week')) = 52;
config(contains(tmp, '2-4x/week')) = 52*3;
config(contains(tmp, 'daily')) = 365;
%impute missing data as median
config(isnan(config)) = nanmedian(config);
X = [X,config];

%social activity yes or no
tmp = demos.R_SI_SocialAct_f6160_2;
config = nan(length(tmp), 1);
config(strcmp(tmp, 'Social Activity')) = 1;
config(strcmp(tmp, 'No Social Activity')) = 0;
%impute missing data as mean
config(isnan(config)) = nanmedian(config);
X = [X,config];

%how often can you confide. Convert into number of times per year
tmp = demos.R_SI_AbletoConfide_f2110_2;
config = nan(length(tmp), 1);
config(contains(tmp, 'Never')) = 0;
config(contains(tmp, 'few months')) = 3;
config(contains(tmp, '1x/month')) = 12;
config(contains(tmp, '1x/week')) = 52;
config(contains(tmp, '2-4x/week')) = 52*3;
config(contains(tmp, 'Daily')) = 365;
%impute missing data as median
config(isnan(config)) = nanmedian(config);
X = [X,config];

%switch to negative so that higher values indicate social ISOLATION
X = -X;

%social isolation score within the first principal component
[~, SCORE] = pca(zscore(X));
social_isolation = SCORE(:,1);

%BMI
BMI = demos.R_BMI_f21001_2; BMI(isnan(BMI)) = nanmedian(BMI);

%Alcohol intake
alcohol = demos.R_AlcoholIntakeCompUnitsPerMonth_final_2; alcohol(isnan(alcohol)) = nanmedian(alcohol);



%Lack of secondary education
tmp = demos.R_RF_DemLC_1_LackSecondaryEdu;
config = nan(length(tmp), 1);
config(strcmp(tmp, 'Yes')) = 1;
config(strcmp(tmp, 'No')) = 0;
config(isnan(config)) = nanmedian(config); %impute missing data as median
education = config;

%Hearing loss
tmp = demos.R_RF_DemLC_2_Hearing_2;
config = nan(length(tmp), 1);
config(strcmp(tmp, 'Yes')) = 1;
config(strcmp(tmp, 'No')) = 0;
config(isnan(config)) = nanmedian(config); %impute missing data as median
hearing = config;

%Traumatic brain injury
tmp = demos.R_RF_DemLC_3_HeadInjury_2;
config = nan(length(tmp), 1);
config(strcmp(tmp, 'Yes')) = 1;
config(strcmp(tmp, 'No')) = 0;
config(isnan(config)) = nanmedian(config); %impute missing data as median
TBI = config;

%History of depression
tmp = demos.R_RF_DemLC_10_Depression_2;
config = nan(length(tmp), 1);
config(strcmp(tmp, 'Yes')) = 1;
config(strcmp(tmp, 'No')) = 0;
config(isnan(config)) = nanmedian(config); %impute missing data as median
depression = config;

%Physical inactivity
tmp = demos.R_RF_DemLC_11_PhysicalInactivity;
config = nan(length(tmp), 1);
config(strcmp(tmp, 'Yes')) = 1;
config(strcmp(tmp, 'No')) = 0;
config(isnan(config)) = nanmedian(config); %impute missing data as median
physical_inactivity  = config;

%Air pollution
tmp = demos.R_RF_DemLC_12_AirPollution;
config = nan(length(tmp), 1);
config(strcmp(tmp, 'Yes')) = 1;
config(strcmp(tmp, 'No')) = 0;
config(isnan(config)) = nanmedian(config); %impute missing data as median
pollution  = config;

%Smoking history
tmp = demos.R_RF_DemLC_8_SmokingEver_2;
config = nan(length(tmp), 1);
config(strcmp(tmp, 'No')) = 0;
config(contains(tmp, 'Yes')) = 1;
config(isnan(config)) = nanmedian(config); %impute missing data as median
smoking = config;

%Hypertension
hypertension = demos.HypertensionZComposite_2_N;

%Diabetes
diabetes = demos.DiabetesZComposite_2_N;

%AD Polygenic risk score
PRS = demos.standard_prs_for_alzheimers_disease_ad_f26206_0_0;
PRS(isnan(PRS)) = nanmedian(PRS); %impute missing data as median

%Age
age = demos.R_Age_f21003_2;

%Sex
sex = demos.R_Sex_f31;

%Social deprivation score
Townsend = demos.R_Townsend_f189; Townsend(isnan(Townsend)) = nanmedian(Townsend);

%Binary variable (case vs control)
Y = ~strcmp(demos.R_ML_DiagbySess_C42C240Xf41270f20002_Dementia_2, 'No');


%% Individual regression models
% TBI EXCLUDED AS ONLY 9 POSITIVES

close all

%predictor variables
vars = [PRS, hypertension,  diabetes, BMI, alcohol, smoking, depression, education, pollution,  physical_inactivity, hearing, social_isolation]

%With labels for visualisation
labels = {'PRS','Hypertension', 'Diabetes', 'BMI', 'Alcohol', 'Smoking', 'Depression', 'Less education','Pollution', 'Physical inactivity', 'Hearing loss',  'Social isolation'}

%Observation weights so that cases upweighted
weights = ones(length(Y),1);
weights(Y==1) = 1-sum(Y==1)./numel(Y);
weights(Y==0) = 1-sum(Y==0)./numel(Y);


p = []; H = []; beta = []; X = [];
for ivar = 1:size(vars, 2)

    %One single predictor + 3 covariates of no-interest
    X = zscore([vars(:, ivar), age, sex, Townsend]);

    %Fit linear regression model to predict EC index
    mdl = fitglm(X,EC_index, 'CategoricalVars', [3], 'weights', weights)

    beta(ivar) = mdl.Coefficients.Estimate(2);
    SE(ivar) = mdl.Coefficients.SE(2);
    p(ivar) = mdl.Coefficients.pValue(2);
end

%correct for multiple comparisons using Holm-Bonferroni
[p_sorted,p_corrected,H, labels, beta, SE] = UKB_DCM_dem_holmbonf(p, labels, beta, SE)


grey = [0.65 0.65 0.65]
purple = [0.4 0 0.8]
axiswidth = 3;
scatter_size = 30;
fsz = 14;

h = barwitherr(SE, beta)
set(h(1), 'FaceColor', grey)
set(h(1), 'EdgeColor', 'none')
set(gca, 'FontSize', fsz)
set(gca, 'LineWidth', axiswidth)
box off
axis([0 13 -0.025 0.1])
predictors = labels
xticklabels(predictors)
ylabel('Beta')
set(gcf,'position',[250,250,1000,600])

sig = find(p_corrected<0.05)
hold on

scatter(sig, 0.08,1000,purple, '.')
title('Associations between risk factors and effective connectivity change')


%% mediation analysis
%MEDIATOR (effective connectivity)
M = EC_index;

%X VARIABLE (social isolation)
X= social_isolation;

%Y VARIABLE (case vs. control)
Y = ~strcmp(demos.R_ML_DiagbySess_C42C240Xf41270f20002_Dementia_2, 'No');


%COVARIATES
covs = [age, sex, Townsend];

%Observation weights so that cases upweighted
weights = ones(length(Y),1);
weights(Y==1) = 1-sum(Y==1)./numel(Y);
weights(Y==0) = 1-sum(Y==0)./numel(Y);

%Number of permutations
nperm = 1000;

mediation_output = UKB_DCM_dem_mediation(X,M,Y,covs,weights, nperm)


end