% life_test
feOpenLocalCluster;

%%
cd ~/Downloads/data
addpath(genpath(pwd))

%% Build the file names for the diffusion data, the anatomical MRI.
dwiFile       = fullfile(lifeDemoDataPath('diffusion'),'life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz');
dwiFileRepeat = fullfile(lifeDemoDataPath('diffusion'),'life_demo_scan2_subject1_b2000_150dirs_stanford.nii.gz');
t1File        = fullfile(lifeDemoDataPath('anatomy'),  'life_demo_anatomy_t1w_stanford.nii.gz');

%% 
fgFileName    = fullfile(lifeDemoDataPath('tractography'), ...
                'life_demo_mrtrix_csd_lmax10_probabilistic.mat');

fg = fgRead(fgFileName);
small_fg = fgCreate('name', ['small_' fg.name], 'fibers', fg.fibers(1:10));
small_fg.pathwayInfo = fg.pathwayInfo(1:10);
fgWrite(small_fg, fullfile(lifeDemoDataPath('diffusion'), 'small_fg.mat'));

%% 
feFileName    = 'life_build_model_demo_CSD_PROB_small';

%% Build a model of the connectome.
fe = feConnectomeInit(dwiFile, small_fg, feFileName,fullfile(fileparts(fgFileName)),dwiFileRepeat,t1File);
%fe = feConnectomeInit(dwiFile, fgFileName, feFileName,fullfile(fileparts(fgFileName)),dwiFileRepeat,t1File);

%% Fit the model with global weights.
fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));

%% 
rmse   = feGet(fe,'vox rmse');
dwi = dwiLoad(feGet(fe,'dwifile'));
%% make a predicted diffusion data
outputname = 'rRMSEmap_life_build_model_demo_CSD_PROB_small.nii.gz';
    
% retrieve argument form fe structure
coords  = feGet(fe,'roi coords');
xform   = feGet(fe,'xform img 2 acpc');
mapsize = feGet(fe, 'map size');
fiber_w = feGet(fe,'fiber weights');
%% Get rmseR
% rmseM   = feGetRep(fe, 'vox rmse');
%rmseD   = feGetRep(fe, 'vox rmse data');
rmseR   = feGetRep(fe, 'vox rmse ratio');
%% Other numbers
% The woxels returned by a fit of LiFE by voxel/fiber
% w = feGet(fe,'fiberweightsvoxelwise')
%---------
% Predict the diffusion signal for the fiber component
% with the voxel-wise fit of LiFE
% pSig = feGet(fe,'psigfvoxelwise')
% pSig = feGet(fe,'psigfvoxelwise',coords)
% pSig = feGet(fe,'psigfvoxelwise',voxelIndices)
%------
% Predicted signal of fiber alone (demeaned).
pSig = feGet(fe,'pSig fiber');
%---------
% Predict the diffusion signal for the fiber component
% with the voxel-wise fit of LiFE, return the an array of pSigXnVoxel
% pSig = feGet(fe,'psigfvoxelwisebyvoxel')
% pSig = feGet(fe,'psigfvoxelwisebyvoxel',coords)
% pSig = feGet(fe,'psigfvoxelwisebyvoxel',voxelIndices)
%---------
% The fiber and isotropic weights as a long vector
% w = feGet(fe,'fullweights')
%---------
% Return the global R2 (fraction of variance explained) of the full life
% model.
% R2 = feGet(fe,'total r2');
% Residual signal: (fiber prediction - measured_demeaned).
% res = feGet(fe,'res sig fiber');

% %% RMSE off the model
% maxRmse = 90;
% % rmseImg = feReplaceImageValues(nan(mapsize),rmseR,coords);
% rmseImg = feReplaceImageValues(zeros(mapsize),rmseR,coords);
% 
% rmseImg(rmseImg > maxRmse) = maxRmse;
% ni  = niftiCreate('data',rmseImg, 'fname', outputname, ...
%     'qto_xyz',xform, ...
%     'fname','FDM', ...
%     'data_type',class(rmseImg));

% %% the total number of fibers for all the voxels
% % total number
% nFibers = feGet(fe,'totfnum');
% size(nFibers')
% nFibersImg = feReplaceImageValues(zeros(mapsize),nFibers',coords);
% outputname = 'TotalNumFibers_life_build_model_demo_CSD_PROB_small.nii.gz';
% 
% ni  = niftiCreate('data',nFibersImg, 'fname', outputname, ...
%     'qto_xyz',xform, ...
%     'fname','nFibers', ...
%     'data_type',class(nFibersImg));

%% Measured signal in VOI, demeaned, this is the signal used for the fiber-portion of the M model.
% dSig = feGet(fe,'dsigdemeaned',coords');
% 
% %---------
% % Demeaned diffusion signal in each voxel.
% dSigByVoxel = feGet(fe,'dsigdemeaned by Voxel');
% %  dSigByVoxel = feGet(fe,'dsigdemeaned by Voxel',coords');
% % dSigByVoxel = feGet(fe,'dsigdemeaned by Voxel',vxIndex);
% size(dSigByVoxel)

%% Predicted signal
%---------
% Predicted signal by the full model in a set of voxeles.
pSigByVoxel = feGet(fe, 'psigfullvox');
% pSigByVoxel2 = feGet(fe,'psigfullvox',coords');
% pSigByVoxel = feGet(fe, 'pSig full by voxel',voxelIndex);
size(pSigByVoxel)

%%
% Reshape the predicted signal by voxles
    predicted = reshape(feGet(fe,'pSig fiber test'), feGet(fe,'nBvecs'), feGet(fe,'nVoxels'));



%%
pSigByVoxelImg = feReplaceImageValues(zeros(mapsize),pSigByVoxel,coords);
outputname = 'pSigByVoxel_life_build_model_demo_CSD_PROB_small.nii.gz';

ni  = niftiCreate('data',pSigByVoxelImg, 'fname', outputname, ...
    'qto_xyz',xform, ...
    'fname','nFibers', ...
    'data_type',class(pSigByVoxelImg));


%% check generated nifti file
showMontage(ni.data)
colormap 'jet'

%% create dwi 
ni  = niftiCreate('data',nFibersImg, 'fname', outputname, ...
    'qto_xyz',xform, ...
    'fname','nFibers', ...
    'data_type',class(nFibersImg));

% Load original dwi
dwi = dwiLoad(feGet(fe,'dwifile'));

% Example:
ni = dwiCreate('nifti',[],'bvecs',dwi.bvecs,'bvals',dwi.bvals);


%%  
if save_flag,
niftiWrite(ni,outputname);end

return