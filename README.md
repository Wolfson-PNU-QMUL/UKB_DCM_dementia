# UKB_DCM_dementia
This repository provides all the code that was used for our project of fitting spectral dynamic causal models (DCM) to dementia cases' and matched controls' resting-state fMRI (rs-fMRI) data from the UK biobank.

You can use this repository to reproduce the results and figures from our original analysis.

Or you can use this repository if you're planning on applying spectral DCM to a new dataset and/or a different brain network. The code should be fairly easy to modify and re-apply for your own project. 

If you're new to DCM, here are some essential resources to get your started:

https://en.wikibooks.org/wiki/SPM/Parametric_Empirical_Bayes_(PEB)
https://www.sciencedirect.com/science/article/pii/S1053811919305221
https://www.sciencedirect.com/science/article/pii/S1053811919305233

## What is included in this repository?
### Code
Code for the entire analysis pipeline, which takes a raw T1-weighted structural MRI image and raw rs-fMRI data and outputs DCM results and performance of predictive models designed to predict dementia incidence and time until future diagnosis. 
### Data
There are also some data files that include the end-results of our analysis pipeline:

Specifically, the BMA (bayesian model average) results which summarise the group-level effective connectivity differences between dementia cases and controls, and the BMA results which summarise the group level effective connectivity differences, within dementia cases, that vary as a function of time until diagnosis. 

There are also files which include the results of K-fold cross-validated logistic regression to predict dementia incidence from effective connectivity, functional connectivity and volumetric data respectively (EC_classifier_Kfold.mat, FunctConx_classifier_Kfold.mat, volumetric_classifier_Kfold.mat), and files which include the results of K-fold cross-validated linear regression to predict time until dementia diagnosis from effective, connectivity, functional connectivity and volumetric data respectively (EC_prognosticator_Kfold.mat, FunctConx_prognosticator_Kfold.mat, volumetric_prognosticator_Kfold.mat).

There are also two files called LOOCV_EC_classifier_weighted_PreClinOnly_1.mat and LOOCV_EC_classifier_weighted_PreClinOnly_0.mat. These are the outputs of leave-one-out cross-validated logistic regression models trained on effective connectivity to predict dementia incidence. The first file (appended preclin_1) excludes cases with prevalent dementia at the time of data acquisition. The second file (appended preclin_0) includes all subjects. These files are used to derive the effective connectivity (EC) index, which is a subject-specific variable describing how "dementia-like" an individual's DMN effective connectivity pattern is. It is simply the P(dementia) outputted by the cross-validated regression model. This subject-specific variable is used for the modifiable risk factors analysis.

## Assumptions on data structure
The code makes some assumptions about how your MRI data files are organised. 

In the master script (UKB_DCM_dem_pipeline.m) you will need to define the directory where your structural and fMRI data are located

```
fMRI_datadir = '';
sMRI_datadir = '';
```

For the structural MRI data the code will assume that within your directory there is a sub-directory for each participant and that within that sub-directory there is a sub-directory called T1, within which is your T1 image labelled T1.nii. E.g.:

```
[sMRI_datadir '/001/T1/T1.nii']
```
where 001 is the subject ID 

For the fMRI data the code will assume that within your directory there is a sub-directory for each participant and that within that sub-directory there is a sub-directory called fMRI, within which your raw data is labelled rfMRI.nii and a single-band reference image is labelled rfMRI.SBREF.nii E.g.: 

```
[fMRI_datadir '/001/fMRI/rfMRI.nii']
[fMRI_datadir '/001/fMRI/rfMRI_SBREF.nii']
```

If you are a QMUL researcher with access to Apocrita then you can use the direct paths to the stored UKB data as follows:

```
fMRI_datadir = '/data/Wolfson-PNU-dementia/UKB/Imaging_rsfMRI_UKB_DCM/rfMRI_CCDem_2';
sMRI_datadir = '/data/Wolfson-PNU-dementia/UKB/Imaging_rsfMRI_UKB_DCM/sMRI_CCDem_2';
```

The code makes reference to 2 csv files containing raw demographic and brain volumetric data from the UK Biobank which cannot be shared on Github. QMUL researchers with access to Apocrita can find the relevant .csv files at the following paths:

```
/data/Wolfson-PNU-dementia/UKB/Imaging_rsfMRI_UKB_DCM/phenotype_CCDem_2/UKBW_rsfmri_F_DemIncident_F10_280623.csv
/data/Wolfson-PNU-dementia/UKB/Imaging_rsfMRI_UKB_DCM/phenotype_CCDem_2/UKBW_rsfmri_F_DemPrevalent_F10_280623.csv
```

## How to use the code
Start off by looking at the file: UKB_DCM_dem_pipeline.m

This script can execute each step of the analyis pipeline in discrete chunks (functions).

The chunks are as follows:

1) Structural T1 image processing  
UKB_DCM_dem_structural_MRI_preprocess.m pre-processes and skull-strips the T1-weighted structural image.

2) fMRI pre-processing  
UKB_DCM_dem_functional_MRI_preprocess.m pre-processes the fMRI data (realignment, co-registration, normalization and smoothing).

3) Timeseries extraction  
UKB_DCM_dem_extract_timeseries.m extracts BOLD timeseries from the brain regions that will be used in the subsequent DCM analysis. 
This is a step where you will need to edit the code if you're interested in a different brain network. You can actually define the labels and MNI co-ordinates you're interested in by editing the function UKB_DCM_dem_ROI_specify.m

4) Fit a fully connected DCM to each individual subject  
UKB_DCM_dem_firstlevelDCM.m

5) Identify differences in DCM parameters between cases and controls and classify
UKB_DCM_dem_EC_classifier.m does two things. Firstly it fits a group level PEB (parametric empirical Bayes) model to all of the individual DCMs (cases and controls) and then prunes the group level PEB using BMA (Bayesian model averaging). Secondly it uses the resultant parameters that are different between cases and controls to try and predict dementia incidence using leave-one-out cross-validated elastic-net logistic regression.

6) Identify differences in DCM parameters as a function of time until diagnosis and predict
UKB_DCM_dem_EC_prognosticator.m does two things. Firstly it fits a group level PEB model to the DCMs of the dementia cases and then prunes the group level PEB using BMA. Secondly it uses the resultant parameters that vary as a function of time until diagnosis to try and predict time until diagnosis using leave-one-out cross-validated elastic-net linear regression.

7) Analysis of modifiable risk factors  
UKB_DCM_dem_ModRiskFactors.m Extracts data on modifiable risk factors and runs a number of independent regression analyses to search for associations between risk factors and "EC index" (EC index is simply the probability of dementia outputted by the classifier in step 5). This function also runs a mediation analysis to test whether EC index mediates the association between social isolation and dementia incidence.

8) Volumetric analysis
UKB_DCM_dem_volumetric.m uses regularised logistic regression to predict dementia incidence using structural MRI data and then uses regularised linear regression to predict time to dementia diagnosis using structural MRI.

9) Functional connectivity analysis
UKB_DCM_dem_functional_connectivity.m uses regularised logistic regression to predict dementia incidence using DMN functional connectivity features and then uses regularised linear regression to predict time to dementia diagnosis using DMN functional connectivity features.



## How to generate figures
Here is some code that you can use to visualise your DCM results using MATLAB

First we can try the simple thing of making a connectivity matrix

Firstly load any BMA (Bayesian Model Average) file. This is the output of the SPM function spm_dcm_peb_bmc.m

```load('BMA_EC_classifier.mat')```

See how many ROIs/nodes were in this DCM

```nR = sqrt(numel(BMA.Pnames));```

Work out which elements are of interest. In this case, the BMA has pruned from an original total of 500 parameters. There were 100 connectivity paramers (10 ROIs fully connected to each other). Then there were 5 columns in the between-subjects design matrix. Therefore a total of 500 parameters. Typically, the design matrix is designed such that the SECOND column is the one we're interested in. The subsequent columns are covariates of no-interest. The first column usually represents a group mean or baseline. In this case, the second column represents differences in connectivity parameters between cases and controls. We are therefore interested in parameters 101 to 200. Define the start and end element of this subset of parameters.

```
start = (nR^2)+1;
finish = start-1+(nR^2);
parameters = BMA.Ep(start:finish) %The parameters values are in a field labelled .Ep
```

Reshape the data into a nR x nR matrix (in this case it's 10 x 10)

```
parameters = reshape(parameters, nR,nR)
```

Identify which parameters have a high posterior probability (99%) of being non-zero

```
Pp = reshape(BMA.Pp(start:finish),nR,nR)>0.99; %Stored in a field labelled .Pp
```

Exclude parameters which don't have a very strong evidence of being non-zero
```
X = parameters.*Pp;
```

Label your ROIs for visualisation (the order here is very important. It must be the same order you used when you first fit your DCM)

```
ROIs = {'PRC', 'amPFC', 'vmPFC', 'dmPFC', 'lIPC', 'rIPC', 'lLTC', 'rLTC', 'lPHF', 'rPHF'}';
```

Use the following function to plot the connectivity matrix
```
PlotConnectMatrix(X, ROIs)
```
The connectivity matrix will look something like this:

![ConnectivityMatrixExample](https://github.com/Wolfson-PNU-QMUL/UKB_DCM_dementia/assets/142310228/4fe250a9-9371-4349-8d1c-d395e8072eab)

Connectivity matrices are clear and concise but not particularly anatomically intuitive. 

An alternative is plot the connectivity parameters in anatomical space within a glass brain.

Here is some code that will let you do that. 

To avoid the images being too crowded the code will make the following different brain-plots. 

1) Strengthened excitatory connections
2) Strengthened inhibitory connections
3) Attenuated excitatory connections
4) Attenuated inhibitory connections

Of course you could edit this code to plot all of the connectivity changes in one single glass brain

Firstly, as before, load any BMA (Bayesian Model Average) file (this is the output of the SPM function spm_dcm_peb_bmc.m

```load('BMA_EC_classifier.mat')```

Now let's initialise some graphical settings (which you can play around with)

```
sphere_size = 8.5 
brain_alpha = 0.05 
isovalue = 25
braincol = [0.8,0.5,1]
```

Select whether you want to view alterations of excitatory or inhibitory connections

```
valence = 1 %excitatory. Set to 2 to switch it to inhibitory
```

Select whether you want to view attenuated connections or strengthened connections

```
direction = 1 %attenuated connections. Set to 2 to switch it to strengthened connections
```

Identify which connections (in baseline/control condition) are excitatory, inhibitory or null

```
nR = sqrt(numel(BMA.Pnames));
Ep_ctrl = sign(full(BMA.Ep(1:nR^2)).*(full(BMA.Pp(1:nR^2))>0.99));
```

Let's define the parameters of interest (connectivity changes) using the same code we used above

```
start = (nR^2)+1;
finish = start-1+(nR^2);
Ep_delta = full(BMA.Ep(start:finish)).*(full(BMA.Pp(start:finish))>0.99);
```

Now we need to select which of these conncetivity changes we are going to plot in this specific glass brain

```
if valence ==1 %modifications to excitatory connections
    if direction ==1
        inds = find(Ep_ctrl>0 & Ep_delta<0); %attenuated excitatory connections
    elseif direction ==2
        inds = find(Ep_ctrl>=0 & Ep_delta>0); %strengthened excitatory connections
    end
elseif valence == 2 %modifications to inhibitory connections
    if direction ==1
        inds = find(Ep_ctrl<0 & Ep_delta>0); %attenuated inhibitory connections
    elseif direction ==2
        inds = find(Ep_ctrl<=0 & Ep_delta<0); %strengthened inhibitory connections
    end
end

% we will also want the parameter values of all the other connections that we're not plotting later on
all_parm_vals = Ep_delta(find(Ep_delta)); 
```

Set non-relevant elements to 0. Reshape it into a matrix.
```
Ep_delta(setdiff([1:numel(Ep_delta)], inds)) = 0;
dat = reshape(Ep_delta,10,10);
```

Get the ROI labels along with MNI co-ordinates. For this project using the DMN we can get those using this function:

```
[ROI_list] = UKB_DCM_dem_ROI_specify();
xyz = cell2mat(ROI_list(:,2)); %MNI co-ordinates for each ROI
```

Define what colour you want each ROI to be using RGB coding. For example: 

```
grey = [0.4, 0.4, 0.4];
red = [0.7,0.1,0.1];
blue = [0.1 0.1 0.7];
yellow = [0.4 0.4 0.2];
purple = [0.6 0.1 0.8];
orange = [0.4 0.2 0];
green = [0.2 0.5 0.2];
facecols = {grey,red,blue,yellow,purple,purple,orange,orange,green,green}; %needs to be the same length as the number ROIs
```

Now plot the relevant ROIs as spheres 

```
retained_nodes = find(sum(dat,1)+sum(dat,2)');
figure
for ixyz = retained_nodes
    [x_range,y,z] = sphere(100);
    x_range = sphere_size*x_range+(xyz(ixyz,1));
    y = sphere_size*y+(xyz(ixyz,2));
    z = sphere_size*z+(xyz(ixyz,3));
    surf(x_range,y,z, 'FaceAlpha', 1, 'EdgeColor', 'none','FaceColor',facecols{ixyz});
    camlight
    axis equal
    hold on
end
```

Now connect the spheres using 3D cylinders 
The width of the cylinder will be proportional to the magnitude of the connectivity change
The colour of cylinder will match the colour of the ROI from which the connection originates
Strengthened connections which be plotted as solid cylinders
Attenuated connections will be plotted as dashed cylinders, using a custom-built function you can download from this repository called dashed_cylinder.m

```
hold on
drawn = []; 
width_spectrum = linspace(sphere_size/5,sphere_size,100); %different line widths
for iparm = 1:length(inds)

    value = Ep_delta(inds(iparm));

    [i,j] = ind2sub([nR,nR],inds(iparm));
    line = [xyz(i,:);xyz(j,:)];
    col = facecols{j};

    % displace line slightly if it connects
    % 2 nodes that are already reverse-connected
    % this is to aid visualisation
    if iparm>1
        if ~isempty(find(squeeze(sum(sum([flipud(line)==drawn],2)) == 6)))
            line = line + 2;
        end
    end

    width_intensity(iparm) = max(round(100.*(abs(value)-min(abs(all_parm_vals)))./max(abs(all_parm_vals))),1);
    width(iparm) = width_spectrum(width_intensity(iparm));

    %draw a line (cylinder) connecting the nodes
    if direction ==1 %attenuated
        [X,Y,Z] = dashed_cylinder(0.5*width(iparm), 100, line(1,:),line(2,:), 6);
        for i = 1:size(X,3)
            surf(X(:,:,i),Y(:,:,i),Z(:,:,i), 'FaceColor', col, 'EdgeColor', 'none');
            camlight
            axis equal
            hold on
        end
    elseif direction == 2 %strengthened
        [X, Y, Z] = cylinder2P(0.5*width(iparm), 100, line(1,:), line(2,:));
        surf(X,Y,Z, 'FaceColor', col, 'EdgeColor', 'none');
        camlight
        axis equal
        hold on
    end

    drawn = cat(3,drawn,line);

end
```

Finally we overlay a glass brain. 

Load the example brain image which you can download from this repository

```
brain = niftiread(['wSkullstr_biascor_structural.nii']);
```

A bit of fiddling to make sure that the co-ordinates match MNI space

```
x_range = [-(size(brain,1)-1)./2:(size(brain,1)-1)./2]';
y_range = [-111:size(brain,2)-111-1]';
z_range = [-70:size(brain,3)-70-1]';

```

Plot the brain

```
[y, x, z] = meshgrid(y_range, x_range, z_range);

patch(isosurface(x,y,z,brain,isovalue),...
    'FaceAlpha', brain_alpha, ...
    'FaceColor', braincol,...
    'EdgeColor', 'none');
```

Find a good viewing angle. You might want to plot the figure twice, in two separate viewing angles

```
angle = 1 %this will give you a posterior to anterior view. Set to 2 to get an anterior to posterior view

if angle ==1
    view([-45.55 22.57]);
elseif angle == 2
    view([-128.62 29.11]);
end
```


Now just some text/axis formatting

```
xlh = xlabel('x', 'FontWeight','bold');
ylh = ylabel('y', 'FontWeight','bold');
zlabel('z','FontWeight', 'bold', 'Rotation', 180)

if angle ==1
    xlh.Position(2) = xlh.Position(2) + abs(xlh.Position(2) * 0.1);
    xlh.Position(1) = xlh.Position(1) + abs(xlh.Position(1) * 5);
    ylh.Position(2)  = ylh.Position(2) + abs(ylh.Position(2) * 1.5);
    ylh.Position(1)  = ylh.Position(1) + abs(ylh.Position(1) * 0.2);
elseif angle ==2
    xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2) * 0.1);
    xlh.Position(1) = xlh.Position(1) + abs(xlh.Position(1) * 2.5);
    ylh.Position(2)  = ylh.Position(2) - abs(ylh.Position(2) * 20);
    ylh.Position(1)  = ylh.Position(1) + abs(ylh.Position(1) * 0.2);
end

set(gca, 'FontSize', 14); set(gca, 'LineWidth', 3);
```

Here's an example of some attenuated connections, viewed from angles 1 and 2 (as defined above):

<img width="1012" alt="AttenuatedConnectionsExample" src="https://github.com/Wolfson-PNU-QMUL/UKB_DCM_dementia/assets/142310228/0639c9a8-a097-41d9-8247-7fffc7d0d48c">


And here's an example of some strengthened connections:

<img width="1012" alt="StrengthenedConnectionsExample" src="https://github.com/Wolfson-PNU-QMUL/UKB_DCM_dementia/assets/142310228/115ad5f7-5bb2-4148-ba39-1fda6f025900">



