function pNifti = runScript

%% Use Franco's example data

try lifeDemoDataPath
catch err
    error('Add life data to your path');
end

%% Build the file names for the diffusion data, the anatomical MRI.

% Here are the original dMRI data
dwiFile       = fullfile(lifeDemoDataPath('diffusion'),'life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz');
dwiFileRepeat = fullfile(lifeDemoDataPath('diffusion'),'life_demo_scan2_subject1_b2000_150dirs_stanford.nii.gz');
t1File        = fullfile(lifeDemoDataPath('anatomy'),  'life_demo_anatomy_t1w_stanford.nii.gz');

%% Create dwi structure

nifti = niftiRead(dwiFile);
bvecs =   dlmread(fullfile(lifeDemoDataPath('diffusion'),'life_demo_scan1_subject1_b2000_150dirs_stanford.bvecs'));
bvals =   dlmread(fullfile(lifeDemoDataPath('diffusion'),'life_demo_scan1_subject1_b2000_150dirs_stanford.bvals'));

dwi   = dwiCreate('nifti',nifti,'bvecs',bvecs','bvals',bvals');

%% load fg file
fgFileName    = fullfile(lifeDemoDataPath('tractography'), ...
    'life_demo_mrtrix_csd_lmax10_probabilistic.mat');

fg = fgRead(fgFileName);

%% Build a model of the connectome

feFileName    = 'life_build_model_demo_CSD_PROB';
fe = feConnectomeInit(dwiFile, fg, feFileName,fullfile(fileparts(fgFileName)),dwiFileRepeat,t1File);

%% Fit the model with global weights.

fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));

%% Predicted diffusion nifti
pNifti =  predictDWInifti(fe, dwi);
% save the file
niftiWrite(pNifti,pNifti.fname);