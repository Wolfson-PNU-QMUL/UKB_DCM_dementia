% Code written by Sam Ereira (s.ereira@qmul.ac.uk) - last updated 2/3/24
% Adapted from Almgren et al. https://github.com/halmgren/Pipeline_preprint_variability_reliability_DCM_rsfMRI
% This code extracts BOLD timecourse from several pre-defined ROIs

function UKB_DCM_dem_extract_timeseries(subjects, fMRI_datadir, maxFD)

funct_ID = '_20227_2_0'; %same for all subjects

%ROI labels and MNI co-ordinates
[ROI_list] = UKB_DCM_dem_ROI_specify;

spm('Defaults','fMRI')

failed_subjs = []; retained_subjs = [];
for isj = 1:numel(subjects)
    tic
    EID = subjects(isj);

    fprintf(['Extracting ROI timeseries for subject: ' num2str(EID) ' (' num2str(isj) ' of ' num2str(numel(subjects)) ')\n'])
    
    funct_data_path = [fMRI_datadir num2str(EID) funct_ID];
    regressors_directory = [funct_data_path '/fMRI/regressors'];
    SPM_file_directory = [funct_data_path '/fMRI/SPM'];

    % skip if timeseries extraction is already complete for this subject
    if exist([SPM_file_directory filesep 'T_ext_flag.mat'])
        load([SPM_file_directory filesep 'T_ext_flag.mat']);
        if strcmp(T_ext_flag, 'failure')
            failed_subjs = [failed_subjs; EID];
        else
            retained_subjs = [retained_subjs; EID];
        end
        clear T_ext_flag
        fprintf('Skipping subject - already done\n')
        continue
    end

    
    %skip subject if motion was too excessive
    load([regressors_directory '/Framewise_Displacement.mat'], 'FD');
    if max(FD) >= maxFD
        fprintf('Skipping subject - motion too excessive\n')
        failed_subjs = [failed_subjs; EID];
        continue
    end

    %try/catch statement to pick up subjects where preprocessing/timeseries extraction failed
    try 
        matlabbatch=[];
        %add this subject's pre-processed fMRI data to the MATLAB search path
        addpath(funct_data_path);
        pre_proc_rfMRI = [funct_data_path '/fMRI/swr_rfMRI.nii'];

        %extract TR from subject header. Should be 0.735 for UKB
        X = niftiinfo(pre_proc_rfMRI);
        TR = X.raw.pixdim(5); clear X

        if round(TR,3) ~= 0.735
            warning('TR read from header is %s...UKB TR should be 0.735', num2str(round(TR,3)))
        end

        %Create regressors (DCT, motion, CSF, WM)
        [nifti_images, DCT] =UKB_DCM_dem_BuildRegressors(funct_data_path, pre_proc_rfMRI, TR);

        %Create directory
        SPM_file_directory = [funct_data_path '/fMRI/SPM'];
        mkdir(SPM_file_directory);

        matlabbatch{1}.spm.stats.fmri_spec.dir = {SPM_file_directory};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(nifti_images);
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});

        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {
            [regressors_directory '/dct_set.mat'] %discrete cosine basis set
            [regressors_directory '/VOI_CSF_signal_1.mat'] %CSF signal
            [regressors_directory '/VOI_WM_signal_1.mat'] %WM signal
            [regressors_directory '/rp.mat'] %6 realignment (motion) parameters
            };

        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

        %%fMRI model estimation
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

        %%Define contrasts
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'Effects_of_interest';
        matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = eye(size(DCT,2)); %F contrast over whole cosine set
        matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.delete = 0;

        spm_jobman('run',matlabbatch);

        % Extract Timeseries from each ROI
        thresholds=nan(1,size(ROI_list,1));

        for VOI_number=1:length(ROI_list) %loop over ROIs
            clear matlabbatch;
            cd(SPM_file_directory)
            matlabbatch{1}.spm.util.voi.spmmat(1) = {['SPM.mat']};
            matlabbatch{1}.spm.util.voi.adjust = 1;
            matlabbatch{1}.spm.util.voi.session = 1;
            matlabbatch{1}.spm.util.voi.name = ROI_list{VOI_number,1};
            matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat(1) = {['SPM.mat']};
            matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = 1;
            matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
            matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
            matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.001;
            matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
            matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
            matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = ROI_list{VOI_number,2};
            matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = 10; %creates 10 mm sphere centered around pre-defined co-ordinates
            matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
            matlabbatch{1}.spm.util.voi.roi{3}.sphere.centre = ROI_list{VOI_number,2};
            matlabbatch{1}.spm.util.voi.roi{3}.sphere.radius = 8; %creates 8 mm sphere centered around peak of F-contrast within the previous 10 mm sphere
            matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.spm = 1;
            matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.mask = 'i2';

            %voxels only included if they exceed the threshold (i1) AND lie
            %within the fixed 10 mm sphere (i2) AND lie within the mobile 8 mm
            %sphere (i3)
            matlabbatch{1}.spm.util.voi.expression = 'i1&i2&i3';

            count=0;
            err_count=0;
            while count==err_count %Make threshold more liberal if no surviving voxels (up to max of 0.05)

                if count==0
                    matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh;
                elseif count==1
                    cd(SPM_file_directory)
                    matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.01;
                elseif count==2
                    cd(SPM_file_directory)
                    matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.05;
                elseif count>=3
                    matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = nan; %These sessions are excluded
                    break
                end
                spm_jobman('run',matlabbatch);

                if ~exist([SPM_file_directory '/VOI_' ROI_list{VOI_number, 1} '_1.mat'])
                    err_count=err_count+1;
                end

                count=count+1;

            end

            thresholds(VOI_number)=matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh;
        end

        save('thresholds.mat', 'thresholds')

        clear matlabbatch

        rmpath(funct_data_path); %to conserve memory

        fprintf(['COMPLETED timeseries extraction for subject: ' num2str(EID) ' (' num2str(isj) ' of ' num2str(numel(subjects)) ')\n'])
        fprintf(['Run time: ' num2str(toc) ' seconds\n'])

        %save flag to show that timeseries extraction is complete
        T_ext_flag = 'success'; 
        save([SPM_file_directory filesep 'T_ext_flag.mat'], 'T_ext_flag');
    catch
        fprintf(['FAILED timeseries extraction for subject: ' num2str(EID) ' (' num2str(isj) ' of ' num2str(numel(subjects)) ')\n'])
        fprintf(['Run time: ' num2str(toc) ' seconds\n'])

        %save flag to show that timeseries extraction is complete
        T_ext_flag = 'failure';
        save([SPM_file_directory filesep 'T_ext_flag.mat'], 'T_ext_flag');
    end

end

end
