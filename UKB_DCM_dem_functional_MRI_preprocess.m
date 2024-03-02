% Code written by Sam Ereira (s.ereira@qmul.ac.uk) - last updated 2/3/24
% Adapted from Almgren et al. https://github.com/halmgren/Pipeline_preprint_variability_reliability_DCM_rsfMRI
% This code pre-processes the UKB rs-fmri data
% No slice timing correction needed due to short TR (multiband sequence)

function UKB_DCM_dem_functional_MRI_preprocess(subjects, fMRI_datadir, sMRI_datadir, maxFD)

struct_ID = '_20252_2_0'; %same for all subjects
funct_ID = '_20227_2_0'; %same for all subjects

spm('Defaults','fMRI')

for isj = 1:numel(subjects)
    tic
    EID = subjects(isj);
    matlabbatch=[];

    fprintf(['Running functional preprocessing for subject: ' num2str(EID) ' (' num2str(isj) ' of ' num2str(numel(subjects)) ')\n'])

    %add this subject's data to the search path
    struct_data_path = [sMRI_datadir num2str(EID) struct_ID];
    funct_data_path = [fMRI_datadir num2str(EID) funct_ID];
    addpath(struct_data_path);
    addpath(funct_data_path);


    %skip if subject already pre-processed
    if exist([funct_data_path filesep 'fMRI/swr_rfMRI.nii'])>0
        fprintf('Done already - skipping...\n')
        continue
    end

    %skip if subject partially pre-processed and already aborted due to
    %high motion
    regressors_directory = [funct_data_path '/fMRI/regressors'];
    if exist([regressors_directory filesep 'Framewise_Displacement.mat'])
        load([regressors_directory filesep 'Framewise_Displacement.mat']);
        if max(FD)>maxFD
            clear FD
            continue
        end
    end

    %Try realigning using single-band reference scan. 
    %If framewise displacement is >2.4 mm then re-try using mean EPI
    %instead
    % If framewise displacement is still >2.4 mm then abort preprocessing
    % and exclude subject
    realigning =1;
    SBREF = 1;
    while realigning
        %Label data files
        rfMRI_file = [funct_data_path '/fMRI/rfMRI.nii'];
        if SBREF ==1
            functional_reference_scan = [funct_data_path '/fMRI/rfMRI_SBREF.nii,1'];
            %Realign to single-band reference image
            data = {[functional_reference_scan;(cellstr(spm_select('expand',rfMRI_file)))]};
        elseif SBREF == 0
            functional_reference_scan = [funct_data_path '/fMRI/meanrfMRI.nii'];
            %Don't realign to single-band reference image 
            data = {[(cellstr(spm_select('expand',rfMRI_file)))]};
        end

        structural_reference_scan = [struct_data_path '/T1/UKBDCM_PreProc/Skullstr_biascor_structural.nii,1'];



        %%%%%%%%%%%%%%%%%%%%%%%
        %Spatial realignment
        %%%%%%%%%%%%%%%%%%%%%%%
        matlabbatch{1}.spm.spatial.realign.estwrite.data = data;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r_';

        spm_jobman('run',matlabbatch);
        clear matlabbatch

        if SBREF == 1
            rp=dlmread([funct_data_path '/fMRI/rp_rfMRI_SBREF.txt']);

            %subtract intial values as this is the reference image (not chronologically first volume)
            rp = rp- repmat(rp(2, :), size(rp,1), 1); 
            R = rp(2:end, :);
        elseif SBREF == 0 
            %if just used mean image instead of SBREF then don't need to subtract
            rp=dlmread([funct_data_path '/fMRI/rp_rfMRI.txt']); R = rp;
        end

        mkdir(regressors_directory);
        save([regressors_directory '/rp.mat'], 'R');

        %Calculate and save framewise displacement
        radius=50; rp(:,4:6)=rp(:,4:6)*radius; rp_diff=diff(rp);
        %add zero as first element (see Power et al., 2014)
        rp_zero=[zeros(1,size(rp_diff,2)); rp_diff];
        FD=sum(abs(rp_zero),2);
        save([regressors_directory '/Framewise_Displacement.mat'], 'FD');

        if max(FD) < maxFD
            fprintf('Framewise displacement is satisfactory\n')
            realigning =0; %exit while loop
        elseif max(FD)>=maxFD && SBREF == 1
            fprintf('Framewise displacement is too high, reattempting without SBREF\n')
            SBREF =0; %stay in while loop and redo realignment
        elseif max(FD)>=maxFD && SBREF == 0
            fprintf('Framewise displacement is still too high, aborting this subject...\n')
            realigning = 0; %exit while loop
        end
    end

    if max(FD)>=maxFD
        continue %skip this subject if FD too high
    end
 
    %%%%%%%%%%%%%%%%%%%%%%%
    %Coregistration (to first anatomical scan)
    %%%%%%%%%%%%%%%%%%%%%%%
    r_rfMRI_file = [funct_data_path '/fMRI/r_rfMRI.nii'];
    realigned_volumes = cellstr(spm_select('expand',r_rfMRI_file));

    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {structural_reference_scan};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {functional_reference_scan};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = realigned_volumes;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    spm_jobman('run',matlabbatch);
    clear matlabbatch

    %%%%%%%%%%%%%%%%%%%%%%%
    %Normalize
    %%%%%%%%%%%%%%%%%%%%%%%

    deformation_field = {[struct_data_path '/T1/UKBDCM_PreProc/y_T1.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = deformation_field;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = realigned_volumes;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
        78 76 85];

    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2.4, 2.4, 2.4];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

    spm_jobman('run',matlabbatch);
    clear matlabbatch

    %%%%%%%%%%%%%%%%%%%%%%%
    %Smooth
    %%%%%%%%%%%%%%%%%%%%%%%

    wr_rfMRI_file = [funct_data_path '/fMRI/wr_rfMRI.nii'];
    normalised_volumes = cellstr(spm_select('expand',wr_rfMRI_file));

    matlabbatch{1}.spm.spatial.smooth.data = normalised_volumes;
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';

    spm_jobman('run',matlabbatch);
    clear matlabbatch

    rmpath(funct_data_path); %to conserve memory
    rmpath(struct_data_path); %to conserve memory

    %save a flag file which indicates whether or not SBREF was used for
    %realignment
    save([funct_data_path filesep 'SBREF_flag.mat'],'SBREF')

    fprintf(['COMPLETED functional preprocessing for subject: ' num2str(EID) ' (' num2str(isj) ' of ' num2str(numel(subjects)) ')\n'])
    fprintf(['Run time: ' num2str(toc) ' seconds\n'])
end

end