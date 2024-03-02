% Code written by Sam Ereira (s.ereira@qmul.ac.uk) - last updated 2/3/24
% Adapted from Almgren et al. https://github.com/halmgren/Pipeline_preprint_variability_reliability_DCM_rsfMRI
% This code fits a fully connected DCM to pre-extracted BOLD timeseries from
% DMN nodes from UK biobank data
% Does a quality check first and only fits DCM if motion is OK and there
% are actually suprathreshold voxels within each of the ROIs.


function UKB_DCM_dem_firstlevelDCM(subjects, fMRI_datadir, maxFD)

spm('Defaults','fMRI')

funct_ID = '_20227_2_0'; %same for all subjects

for isj = 1:numel(subjects)
    tic
    EID = subjects(isj);
    fprintf(['Fitting DCM to subject: ' num2str(EID) ' (' num2str(isj) ' of ' num2str(numel(subjects)) ')\n'])
    funct_data_path = [fMRI_datadir num2str(EID) funct_ID];
    regressors_directory = [funct_data_path '/fMRI/regressors'];

    %Skip this subject if DCM already estimated
    if exist([funct_data_path '/fMRI/Full_DCM/Diagnostics.mat'])>0
        rmpath(funct_data_path) % to conserve memory
        fprintf('DCM already done - skipping subject \n')
        continue;
    end

    %skip subject if motion was too excessive
    load([regressors_directory '/Framewise_Displacement.mat'], 'FD');
    if max(FD) >= maxFD
        fprintf('Skipping subject - motion too excessive\n')
        clear FD
        continue
    end

    %Check if surviving voxels in all ROIs
    ROI_list = UKB_DCM_dem_ROI_specify; nROI = length(ROI_list);
    ROI_sizes = nan(length(ROI_list),1);
    VOI_data_path = [fMRI_datadir num2str(EID) funct_ID '/fMRI/SPM'];
    addpath(VOI_data_path);
    for iROI = 1:length(ROI_list)
        ROI = ROI_list(iROI,1);
        if exist(['VOI_' ROI{1} '_1.mat'])
            load(['VOI_' ROI{1} '_1.mat'])
            ROI_sizes(iROI) = size(xY.XYZmm,2);
        end
    end
    rmpath(VOI_data_path)

    %must be at least 1 suprathreshold voxel in all ROIs
    if sum(ROI_sizes > 0) ~= nROI
        fprintf('No suprathreshold voxels in 1 or more ROIs. Excluding this subject from DCM analysis \n')
        continue
    end

    % Now proceed with DCM fitting if survived quality check


    addpath(funct_data_path);
    pre_proc_rfMRI = [funct_data_path '/fMRI/swr_rfMRI.nii'];
    raw_rfMRI = [funct_data_path '/fMRI/rfMRI.nii'];
    nifti_images = cellstr(spm_select('expand',pre_proc_rfMRI));
    length_scan=size(nifti_images,1);

    %extract TR and TE (in seconds) from subject header. 
    %Should be 0.735 and 0.039 for UKB
    X = niftiinfo(raw_rfMRI);
    TR = double(X.raw.pixdim(5)); 
    TE = str2double(X.raw.descrip(4:5))/1000; clear X

    if round(TR,3) ~= 0.735
        warning('TR read from header is %s...UKB TR should be 0.735', num2str(round(TR,3)))
    end

    if round(TE,3) ~= 0.039
        warning('TE read from header is %s...UKB TE should be 0.039', num2str(round(TE,3)))
    end

    cd([funct_data_path filesep 'fMRI/SPM'])
    tmp = 0;
    for VOI_number=1:size(ROI_list,1)
        ntwrk=ROI_list{VOI_number,1}(1:3);
        node = ROI_list{VOI_number,1}(5:end);
        ntwrk_name{VOI_number,1}=spm_select('FPList',pwd,['^VOI_' ntwrk '_' node '.*\.mat$']);
    end


    DCM = [];

    number_regions=size(ROI_list,1);

    DCM.a=ones(number_regions,number_regions);
    DCM.b=zeros(number_regions,number_regions);
    DCM.c=zeros(number_regions,1);
    DCM.d=zeros(number_regions,number_regions,0);

    DCM.U.u=zeros(length_scan,1);
    DCM.U.name={'null'};

    xY=[];
    for region_number=1:number_regions
        K=load(ntwrk_name{region_number},'xY');
        xY = spm_cat_struct(xY,K.xY);
    end

    DCM.xY=xY;

    DCM.Y.dt=TR;
    DCM.Y.X0=xY(1).X0;
    DCM.Y.Q=spm_Ce(ones(1,number_regions)*length_scan);

    for  region_number = 1:number_regions
        DCM.Y.y(:,region_number)  = xY(region_number).u;
        DCM.Y.name{region_number} = xY(region_number).name;
    end

    DCM.v=length_scan;
    DCM.n=number_regions;

    DCM.TE=TE;
    DCM.delays=ones(number_regions,1)*TR/2;

    DCM.options.nonlinear=0;
    DCM.options.two_state=0;
    DCM.options.stochastic=0;
    DCM.options.centre=0;
    DCM.options.induced=1; %Fit DCM to BOLD cross-spectra

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Estimate DCM files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd([funct_data_path '/fMRI'])

    mkdir(['Full_DCM']); cd(['Full_DCM']);
    rmpath(funct_data_path) % to conserve memory

  
    tic
    DCM=spm_dcm_fit(DCM);
    time_to_converge = toc;

    save('DCM_full_estim.mat','DCM');

    delete([sprintf('DCM_%s',date) '.mat']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Diagnostics on quality of fit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Number of estimable parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DCM = DCM{1};
    qE    = spm_vec(DCM.Ep);
    pE    = spm_vec(DCM.M.pE);
    qC    = DCM.Cp;
    pC    = DCM.M.pC;
    k     = rank(full(pC));
    pC    = pinv(pC);

    D  = trace(pC*qC) + (pE - qE)'*pC*(pE - qE) - spm_logdet(qC*pC) - k;
    N_est_par  = D/log(DCM.v);

    clear qE pE qC pC k D;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %maximum extrinsic connection strength
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Max_conn=max(max(abs(DCM.Ep.A - diag(diag(DCM.Ep.A)))));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Proportion explained variance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PSS   = sum(sum(sum(abs(DCM.Hc).^2)));
    RSS   = sum(sum(sum(abs(DCM.Rc).^2)));
    Expl_var = 100*PSS/(PSS + RSS);

    clear RSS PSS;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Free energy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Free_energy=DCM.F;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Error bars and MAP estimates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    qC  = DCM.Cp;
    qE  = spm_vec(DCM.Ep);
    i=spm_fieldindices(DCM.Ep,'A');
    E=qE(i);
    C=qC(i,i);

    j = 1:size(E,1);

    ci    = spm_invNcdf(1 - 0.05);
    C = diag(C);
    c = ci*sqrt(C(j,:));

    E=spm_unvec(E,zeros(size(DCM.Ep.A,1),size(DCM.Ep.A,1)));
    c=spm_unvec(c,zeros(size(DCM.Ep.A,1),size(DCM.Ep.A,1)));

    Posterior_estimates=E;
    Length_error_bar=c;
    Lower_bound_ci=E-c;
    Higher_bound_ci=E+c;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Satisfactory model convergence?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Model must explain > 10% variance, there must be at least one estimable
    % parameter and at least one extrinsic connection of >1/8 Hz
    converged = Expl_var > 10 && Max_conn > 1/8 && N_est_par >= 1; %1 if converged. 0 if not

    clear qC qE i E C j ci C c;

    save(['Diagnostics.mat'],'time_to_converge' ,'converged', 'N_est_par','Max_conn','Expl_var','Posterior_estimates','Length_error_bar','Lower_bound_ci','Higher_bound_ci','Free_energy');
    clear time_to_converge converged N_est_par Max_conn Expl_var Posterior_estimates Length_error_bar Lower_bound_ci Higher_bound_ci Free_energy;

    fprintf(['Run time: ' num2str(toc) ' seconds\n'])
end

end