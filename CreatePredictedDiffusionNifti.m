function CreatePredictedDiffusionNifti
%% Use Franco's example data

try lifeDemoDataPath
catch
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

dwi   = dwiCreate('nifti',nifti,'bvecs',bvecs,'bvals',bvals);

%% Take the first 10 fibers from the fiber group

% This is intended to run quickly
fgFileName    = fullfile(lifeDemoDataPath('tractography'), ...
    'life_demo_mrtrix_csd_lmax10_probabilistic.mat');

fg = fgRead(fgFileName);
small_fg = fgCreate('name', ['small_' fg.name], 'fibers', fg.fibers(1:10));
small_fg.pathwayInfo = fg.pathwayInfo(1:10);
% fgWrite(small_fg, fullfile(lifeDemoDataPath('diffusion'), 'small_fg.mat'));

clear fg

% small_fg.fibers{ii} are the coordinates in ACPC space

%% Build a model of the connectome.
%
% Running with the issue2 branch on LiFE

feFileName    = 'life_build_model_demo_CSD_PROB_small';
fe = feConnectomeInit(dwiFile, small_fg, feFileName,fullfile(fileparts(fgFileName)),dwiFileRepeat,t1File);
%fe = feConnectomeInit(dwiFile, fgFileName, feFileName,fullfile(fileparts(fgFileName)),dwiFileRepeat,t1File);

%% Fit the model with global weights.

fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));
%% the weights from the LiFE solution

wgts = feGet(fe,'fiber weights');


%% Make an empty nifti file the same size as the original
% pNifti = niftiCreate;
pData = zeros(size(nifti.data));

% duplicate original nifti structure
pNifti      = nifti;
pNifti.data = pData;

% strip .extension
[p,f] = fileparts(pNifti.fname);
[~,f] = fileparts(f);

% give name to pNifti
niftiSet(pNifti,'filepath',[p,f,'_Predicted.nii.gz']);
%% Compute the diffusion tensors for each node in each fiber

% The direction at the node is the average direction between this node and
% each of its neighbors% Isotropic portion of M matrix. Returns the matrix for the full model or

% The diagonal parameters are for a perfect stick, [1 0 0]
Q = feComputeCanonicalDiffusion(fe.fg.fibers, [1 0 0]);  % Q = feGet(fe,''fibers tensors');

%% Add diffusion signal for each fiber coordinate
% We are not sure about which coordinate is the xyz
% We are not sure how to get the S0 value out of the b=0 (non-diffusion
% weighted) image

%%
fe_bvecs2             = bvecs'; %feGet(fe,'bvecs');
fe_bvals2             = bvals'./1000; %feGet(fe,'bvals');

%% nFibers = length(fe.fg.fibers);
for ii = 1:length(fe.fg.fibers);
    
    % These coordinates are within the first fiber, and thus the connectome
    % This is an N x 3 matrix of integers, needed by feGet.
    % This way of getting the fiber coordinates matches the coordinates in the
    % ROI.  But it is obviously ridiculous and if we have to have a difference
    % in the coordinate frames, then there has to be a roi2fg() coordinate
    % transform.  We can't write code like this everywhere.
    fCoords = ((uint32(ceil(fe.fg.fibers{ii})))+1)';
    
    % take S0 image
    S0 = dwiGet(dwi,'b0 image',fCoords); % S0_All = fe.diffusion_S0_img;
    
    % keep indices of the coords in img space.
    indx = sub2ind(dwi.nifti.dim(1:3),fCoords(:,1),fCoords(:,2),fCoords(:,3));
    
    %     sVovel = ind2sub(dwi.nifti.dim(1:3) , indx(jj))
    %% 106 (96 diffusion direction + 10 b0)
    
    for jj = 1:length(indx)
        voxTesor = Q{ii}(jj,:);
        PSig_voxel2{ii,jj} = feComputeSignal(S0(jj), fe_bvecs2, fe_bvals2, voxTesor);
    end
    clear S0
end

%%
pNifti = niftiSet(pNifti,'data',pData);

return

%% -- Just keeping idea



%%
for ii=1:length(nVoxels)
    for jj=1:length(fe_bvecs)
        pData(oneFiber(1,1),oneFiber(2,1),oneFiber(3,1),jj) = ...
            wgts(ii)*feComputeSignal(S0, bvecs', bvals(:), Q{1});
    end
end

pNifti = niftiSet(pNifti,'data',pData);


%% 96 diffusion direction only
fe_bvecs             = feGet(fe,'bvecs');
fe_bvals             = bvals(dwi.bvals ~= 0)'./1000; %feGet(fe,'bvals');

for ii = 1:length(fe.fg.fibers);
    
    % These coordinates are within the first fiber, and thus the connectome
    % This is an N x 3 matrix of integers, needed by feGet.
    % This way of getting the fiber coordinates matches the coordinates in the
    % ROI.  But it is obviously ridiculous and if we have to have a difference
    % in the coordinate frames, then there has to be a roi2fg() coordinate
    % transform.  We can't write code like this everywhere.
    fCoords = ((uint32(ceil(fe.fg.fibers{ii})))+1)';
    
    % take S0 image
    S0 = dwiGet(dwi,'b0 image',fCoords); % S0_All = fe.diffusion_S0_img;
    
    % keep indices of the coords in img space.
    indx = sub2ind(dwi.nifti.dim(1:3),fCoords(:,1),fCoords(:,2),fCoords(:,3));
    
    for jj = 1:length(S0)
        % Once we get the S0 values for this particular voxel, we can compute
        
        voxTesorQ = Q{ii}(jj,:);
        % voxTesorQ = feGet(fe,'voxel tensors',foundVoxels(jj));
        % TensorVoxel = dwiGet(dwi,'tensor image', fCoords(jj,:)); ??
        
        PSig_voxel{ii,jj} = feComputeSignal(S0(jj), fe_bvecs, fe_bvals, voxTesorQ);
    end
    clear S0, clear indx
end

%% Test script
%---------
% Fiber density statistics.
% Computes the fiber density (how many fibers in each voxel)
% We compute the following values:
% (1) The number of fibers in each voxel
% (2) The number of unique fibers with non-zero weights
% (3) The sum of the weights in each voxel
% (4) The mean of the weights in each voxel
% (5) The variance of the weigths in each voxel

FiberStats = feGet(fe,'fiberdensity');

%%
% Keep original index
indx = sub2ind(dwi.nifti.dim(1:3),coords(:,1),coords(:,2),coords(:,3));






