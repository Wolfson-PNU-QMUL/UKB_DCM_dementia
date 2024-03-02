function [ROI_list] = UKB_DCM_dem_ROI_specify()



ROI_list={...

% CORE NETWORK
'DMN_PRC',[2,-58,30],1;...      %% Precuneus/posteriorPCC...ROI from Almgren et al. 2018
'DMN_amPFC',[2,56,-4],2;...           %% mPFC...ROI from Almgren et al. 2018
'DMN_vmPFC',[0 26 18],3;...           %% vmPFC...ROI from Dunn et al. 2014
'DMN_dmPFC',[0 52 26],4;...           %% dmPFC...ROI from Dunn et al. 2014
'DMN_lIPC',[-44,-60,24],5;...         %% inferior parietal...ROI from Almgren et al. 2018
'DMN_rIPC',[54,-62,28],6;...           %%  inferior parietal...ROI from Almgren et al. 2018

%ADDITIONAL NETWORK (TEMPORAL)
'DMN_lLTC', [-60 -24 18],7;...           %% Left lateral temporal cortex (Dunn et al. 2014)
'DMN_rLTC', [60 -24 18],8;...            %% Right lateral temporal cortex (Dunn et al. 2014)
'DMN_lParaHPC', [-28 -40 -12],9;...     %% Left parahippocampal formation (Dunn et al. 2014)
'DMN_rParaPHPC', [28 -40 -12],10;      %% Right parahippocampal formation (Dunn et al. 2014)
};
end