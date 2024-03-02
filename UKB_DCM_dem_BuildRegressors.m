%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Build regressors for DCT, WM, CSF
%Adapted from Almgren et al. https://github.com/halmgren/Pipeline_preprint_variability_reliability_DCM_rsfMRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nifti_images, DCT] =UKB_DCM_dem_BuildRegressors(funct_data_path, pre_proc_rfMRI, TR)

%Construct discrete cosine set (DCT)(adapted from Peter Zeidman (shared via spm archives))

%Extract information regarding number of scans
nifti_images = cellstr(spm_select('expand',pre_proc_rfMRI));
length_scan=size(nifti_images,1); %Usually 490 volumes

%Part adapted from code Peter Zeidman
dct_set = spm_dctmtx(length_scan,length_scan);

upper_limit       = 0.1;
lower_limit       = 1/128;

lower_limit_comp   = fix(2*(length_scan*TR)*lower_limit+1);
dct_set(:,1:lower_limit_comp) = [];

upper_limit_comp   = fix(2*(length_scan*TR)*upper_limit+1);
dct_set(:,upper_limit_comp:end) = [];
R=dct_set;

regressors_directory = [funct_data_path '/fMRI/regressors'];
save([regressors_directory filesep 'dct_set.mat'],'R');

%Make and covert SPM to get compute CSF and WM signal!
matlabbatch{1}.spm.stats.fmri_spec.dir = {regressors_directory};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;

matlabbatch{1}.spm.stats.fmri_spec.sess.scans = nifti_images;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.5;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{3}.spm.util.voi.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.util.voi.adjust = NaN;
matlabbatch{3}.spm.util.voi.session = 1;
matlabbatch{3}.spm.util.voi.name = 'WM_signal';
matlabbatch{3}.spm.util.voi.roi{1}.mask.image(1) = cfg_dep('Model estimation: Analysis Mask', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','mask'));
matlabbatch{3}.spm.util.voi.roi{1}.mask.threshold = 0.5;
matlabbatch{3}.spm.util.voi.roi{2}.sphere.centre = [0 -24 -33];
matlabbatch{3}.spm.util.voi.roi{2}.sphere.radius = 6;
matlabbatch{3}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
matlabbatch{3}.spm.util.voi.expression = 'i1&i2';

matlabbatch{4}.spm.util.voi.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.util.voi.adjust = NaN;
matlabbatch{4}.spm.util.voi.session = 1;
matlabbatch{4}.spm.util.voi.name = 'CSF_signal';
matlabbatch{4}.spm.util.voi.roi{1}.mask.image(1) = cfg_dep('Model estimation: Analysis Mask', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','mask'));
matlabbatch{4}.spm.util.voi.roi{1}.mask.threshold = 0.5;
matlabbatch{4}.spm.util.voi.roi{2}.sphere.centre = [0 -40 -5];
matlabbatch{4}.spm.util.voi.roi{2}.sphere.radius = 5;
matlabbatch{4}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
matlabbatch{4}.spm.util.voi.expression = 'i1&i2';


spm_jobman('run',matlabbatch);
clear matlabbatch

DCT = R;

end
