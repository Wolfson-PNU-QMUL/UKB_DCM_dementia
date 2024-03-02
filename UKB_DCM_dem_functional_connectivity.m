%%% Sam Ereira (s.ereira@qmul.ac.uk) - last updated 2/3/24

function UKB_DCM_dem_functional_connectivity(demos)

subjects = demos.EID;
funct_ID = '_20227_2_0'; %same for all subjects

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% FIRST: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT THE CASE-CONTROL CLASSIFIER %%%%%%%%%%%%%%%%%%%%%%%%
% TO FUNCTIONAL CONNECTIVITY DATA %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute Pearson correlation coefficient between
%timeseries from every pair of ROIs
SimilarityMatrix = []; Y=[];
for isj = 1:length(subjects)
    EID = subjects(isj);
    funct_data_path = ['/data/Wolfson-PNU-dementia/UKB/Imaging_rsfMRI_UKB_DCM/rfMRI_CCDem_2/' num2str(EID) funct_ID];
    load([funct_data_path '/fMRI/Full_DCM/DCM_full_estim']);
    nR = length(DCM{1,1}.xY);
    for roi_i = 1:nR
        for roi_j = 1:nR
            SimilarityMatrix(roi_i,roi_j,isj) = (corr(DCM{1,1}.xY(roi_i).u, DCM{1,1}.xY(roi_j).u));
        end
    end

    %Binary variable (case or control)
    Y(isj,1) =  ~strcmp(demos.R_ML_DiagbySess_C42C240Xf41270f20002_Dementia_2(isj), 'No');
end

%reshape similarity matrix into feature vectors
X=[];
for isj = 1:length(subjects)
    X(isj,:) = reshape(SimilarityMatrix(:,:, isj), 1, nR.^2);
end

%remove self connections (diagonal)
for i_auto = nR:-1:1
    X(:,i_auto+((i_auto-1)*nR)) = [];
end

%Remove symmetrical connections and fisher z transform
X = atanh(unique(round(X',5), 'rows')');

%Fit cross-validated predictive model (elastic-net logistic regression)
[metrics] = UKB_DCM_dem_general_classifier(X,Y,10,5);
dt = string(datetime, 'dd-MM-yy_hh:mm:ss');
save(['FunctConx_classifier_' dt{1} '.mat'], 'metrics')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SECOND: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT THE PROGNOSTICATOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO FUNCTIONAL CONNECTIVITY DATA %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute Pearson correlation coefficient between
%timeseries from every pair of ROIs
SimilarityMatrix = []; Y=[];
sjc=0;
for isj = 1:length(subjects)
    EID = subjects(isj);
    %Only cases
    if ~strcmp(demos.R_ML_DiagbySess_C42C240Xf41270f20002_Dementia_2(find(sum(demos.EID == EID,2))), 'No')
        sjc = sjc+1;

        if EID == 1595913 %%Exclude this subject as unsure when was diagnosed
            continue
        end

        funct_data_path = ['/data/Wolfson-PNU-dementia/UKB/Imaging_rsfMRI_UKB_DCM/rfMRI_CCDem_2/' num2str(EID) funct_ID];
        load([funct_data_path '/fMRI/Full_DCM/DCM_full_estim']);
        nR = length(DCM{1,1}.xY);
        for roi_i = 1:nR
            for roi_j = 1:nR
                SimilarityMatrix(roi_i,roi_j,sjc) = (corr(DCM{1,1}.xY(roi_i).u, DCM{1,1}.xY(roi_j).u));
            end
        end

        %Continuous variable - time to dementia diagnosis
        Y(sjc,1) = str2double(cell2mat(demos.ML_C42C240Xf41270f20002_DementiaAge(isj))) - cell2mat(table2cell(demos(isj,'R_Age_f21003_2')));
    end
end

%reshape similarity matrix into feature vectors
X=[];
for isj = 1:length(Y)
    X(isj,:) = reshape(SimilarityMatrix(:,:, isj), 1, nR.^2);
end

%remove self connections (diagonal)
for i_auto = nR:-1:1
    X(:,i_auto+((i_auto-1)*nR)) = [];
end

%Remove symmetrical connections and fisher z transform
X = atanh(unique(round(X',5), 'rows')');

%Fit cross-validated predictive model (elastic-net logistic regression)
[metrics] = UKB_DCM_dem_general_prognosticator(X,Y,10,5);
dt = string(datetime, 'dd-MM-yy_hh:mm:ss');
save(['FunctConx_prognosticator_' dt{1} '.mat'], 'metrics')

end