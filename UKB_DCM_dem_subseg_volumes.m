function [DAT,names] = UKB_DCM_dem_subseg_volumes(demos)

subjects = demos.EID;
subseg = readtable('/data/Wolfson-PNU-dementia/UKB/Imaging_rsfMRI_UKB_DCM/phenotype_CCDem_2/UKBW_78867r672482_FSSubSegandCognition.csv');
[C,IA,IB] = intersect(subseg.EID59138,subjects, 'stable');
subseg = subseg(IA,:);

%%

DAT = [];
names = {};

tmp = table2array(subseg(:,'volume_of_ca1body_left_hemisphere_f26622_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_ca1body_left'];

tmp = table2array(subseg(:,'volume_of_ca1body_right_hemisphere_f26644_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_ca1body_right'];

tmp = table2array(subseg(:,'volume_of_ca1head_left_hemisphere_f26626_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_ca1head_left'];

tmp = table2array(subseg(:,'volume_of_ca1head_right_hemisphere_f26648_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_ca1head_right'];

tmp = table2array(subseg(:,'volume_of_ca3body_left_hemisphere_f26632_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_ca3body_left'];


tmp = table2array(subseg(:,'volume_of_ca3body_right_hemisphere_f26654_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_ca3body_right'];

tmp = table2array(subseg(:,'volume_of_ca3head_left_hemisphere_f26637_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_ca3head_left'];

tmp = table2array(subseg(:,'volume_of_ca3head_right_hemisphere_f26659_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_ca3head_right'];

tmp = table2array(subseg(:,'volume_of_ca4body_left_hemisphere_f26635_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_ca4body_left'];

tmp = table2array(subseg(:,'volume_of_ca4body_right_hemisphere_f26657_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_ca4body_right'];

tmp = table2array(subseg(:,'volume_of_ca4head_left_hemisphere_f26634_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_ca4head_left'];


tmp = table2array(subseg(:,'volume_of_ca4head_right_hemisphere_f26656_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_ca4head_right'];

tmp = table2array(subseg(:,'volume_of_hippocampaltail_left_hemisphere_f26620_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_hippocampaltail_left'];

tmp = table2array(subseg(:,'volume_of_hippocampaltail_right_hemisphere_f26642_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_hippocampaltail_right'];

tmp = table2array(subseg(:,'volume_of_wholehippocampalbody_left_hemisphere_f26639_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_wholehippocampalbody_left'];


tmp = table2array(subseg(:,'volume_of_wholehippocampalbody_right_hemisphere_f26661_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_wholehippocampalbody_right'];


tmp = table2array(subseg(:,'volume_of_wholehippocampalhead_left_hemisphere_f26640_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_wholehippocampalhead_left'];

tmp = table2array(subseg(:,'volume_of_wholehippocampalhead_right_hemisphere_f26662_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_wholehippocampalhead_right'];


tmp = table2array(subseg(:,'volume_of_wholehippocampus_left_hemisphere_f26641_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_wholehippocampus_left'];



tmp = table2array(subseg(:,'volume_of_wholehippocampus_right_hemisphere_f26663_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_wholehippocampus_right'];


tmp = table2array(subseg(:,'volume_of_hippocampalfissure_left_hemisphere_f26624_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_wholehippocampalfissure_left'];



tmp = table2array(subseg(:,'volume_of_hippocampalfissure_right_hemisphere_f26646_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_wholehippocampalfissure_right'];



tmp = table2array(subseg(:,'volume_of_molecularlayerhpbody_left_hemisphere_f26630_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_molecularlayerhpbody_left'];


tmp = table2array(subseg(:,'volume_of_molecularlayerhpbody_right_hemisphere_f26652_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_molecularlayerhpbody_right'];


tmp = table2array(subseg(:,'volume_of_molecularlayerhphead_right_hemisphere_f26651_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_molecularlayerhphead_right'];


tmp = table2array(subseg(:,'volume_of_molecularlayerhphead_left_hemisphere_f26629_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_molecularlayerhphead_left'];

tmp = table2array(subseg(:,'volume_of_parasubiculum_left_hemisphere_f26628_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_parasubiculum_left'];

tmp = table2array(subseg(:,'volume_of_parasubiculum_right_hemisphere_f26650_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_parasubiculum_right'];

tmp = table2array(subseg(:,'volume_of_presubiculumbody_left_hemisphere_f26627_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_presubiculumbody_left'];

tmp = table2array(subseg(:,'volume_of_presubiculumbody_right_hemisphere_f26649_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_presubiculumbody_right'];

tmp = table2array(subseg(:,'volume_of_presubiculumhead_left_hemisphere_f26625_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_presubiculumhead_left'];

tmp = table2array(subseg(:,'volume_of_presubiculumhead_right_hemisphere_f26647_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_presubiculumhead_right'];

tmp = table2array(subseg(:,'volume_of_subiculumbody_left_hemisphere_f26621_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_subiculumbody_left'];


tmp = table2array(subseg(:,'volume_of_subiculumbody_right_hemisphere_f26643_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_subiculumbody_right'];

tmp = table2array(subseg(:,'volume_of_subiculumhead_left_hemisphere_f26623_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_subiculumhead_left'];

tmp = table2array(subseg(:,'volume_of_subiculumhead_right_hemisphere_f26645_2_0'));
DAT = [DAT,cellfun(@(x) str2double(x), tmp)];
names = [names; 'volume_of_subiculumhead_right'];


tmp = table2array(demos(:, 'volume_of_grey_matter_in_parahippocampal_gyrus_anterior_divisio'))
DAT = [DAT,tmp];
names = [names; 'volume_of_parahippocampalgyrus_anterior_left'];

tmp = table2array(demos(:, 'volume_of_grey_matter_in_parahippocampal_gyrus_anterior_divis_1'))
DAT = [DAT,tmp];
names = [names; 'volume_of_parahippocampalgyrus_anterior_right'];


tmp = table2array(demos(:, 'volume_of_grey_matter_in_parahippocampal_gyrus_posterior_divisi'))
DAT = [DAT,tmp];
names = [names; 'volume_of_parahippocampalgyrus_posterior_left'];


tmp = table2array(demos(:, 'volume_of_grey_matter_in_parahippocampal_gyrus_posterior_divi_1'))
DAT = [DAT,tmp];
names = [names; 'volume_of_parahippocampalgyrus_posterior_right'];

