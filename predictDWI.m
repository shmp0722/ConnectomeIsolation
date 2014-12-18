function CreatePredictedDiffusionNifti
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

%% Take the first 10 fibers from the fiber group

% This is intended to run quickly
fgFileName    = fullfile(lifeDemoDataPath('tractography'), ...
    'life_demo_mrtrix_csd_lmax10_probabilistic.mat');

fg = fgRead(fgFileName);
small_fg = fgCreate('name', ['small_' fg.name], 'fibers', fg.fibers(1:10));
small_fg.pathwayInfo = fg.pathwayInfo(1:10);
% fgWrite(small_fg, fullfile(lifeDemoDataPath('diffusion'), 'small_fg.mat'));

clear fg

% transform coordinates in image space form acpc
small_fg_img = dtiXformFiberCoords(small_fg, inv(dwi.nifti.qto_xyz), 'img');

% img_coords  = horzcat(small_fg_img.fibers{:});
% img_coords2 = horzcat(fe.fg.fibers{:}); These are identical

%% Build a model of the connectome.
%
% Running with the issue2 branch on LiFE
%
% If we don't need fiber weights, we can skip this

feFileName    = 'life_build_model_demo_CSD_PROB_small';
fe = feConnectomeInit(dwiFile, small_fg, feFileName,fullfile(fileparts(fgFileName)),dwiFileRepeat,t1File);

%% Fit the model with global weights.

fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));

%% Make an empty nifti file the same size as the original
% pNifti = niftiCreate;
pData = zeros(size(nifti.data));

% duplicate original nifti structure
pNifti      = nifti;
pNifti.data = pData;

% strip .extension
[p,f] = fileparts(pNifti.fname);
[~,f] = fileparts(f);
newFname = [p,f,'_Predicted.nii.gz'];
% give name to pNifti
pNifti.fname = newFname;

%% Here is the predicted signal from the life fit for the whole connectome

Mfiber = feGet(fe,'M fiber');
Miso   = feGet(fe,'M iso');
wgts   =  feGet(fe,'full weights');
pSig   = [Mfiber,Miso]*wgts;

coords = feGet(fe,'roi coords');
nVoxels = size(coords,1);

% Take the pSig values and put them in a data volume like the nifti data
% volume
roiSig = zeros(size(nifti.data));
for cc=1:size(coords,2)
    roiSig(coords(cc,1),coords(cc,2),coords(cc,3),11:end) = pSig((cc-1)*96 + (1:96));
end

% Or, to compare the predicted and observed signal, get the observed signal
% into a vector like pSig
sig  = niftiGet(nifti,'data');
oSig = zeros(size(pSig));
for cc = 1:nVoxels
    oSig((cc-1)*96 + (1:96)) = sig(coords(cc,1),coords(cc,2),coords(cc,3),11:end);
end

%% Debug, what is wrong with the coordinates if anything ...

mrvNewGraphWin;
plot(oSig(:),pSig(:),'.')


%% Figure out how to zero out some of the fibers
newWgts = wgts;
newWgts(3:8) = 0;

pSig2   = [Mfiber,Miso]*newWgts;




%% Compute the diffusion tensors for each node in each fiber

% The direction at the node is the average direction between this node and
% each of its neighbors% Isotropic portion of M matrix. Returns the matrix for the full model or

% The diagonal parameters are for a perfect stick, [1 0 0]
Q = feComputeCanonicalDiffusion(small_fg_img.fibers, [1 0 0]);  % Q = feGet(fe,''fibers tensors');
% Q = feComputeCanonicalDiffusion_voxel(img_coords, [1 0 0]);
%% Add diffusion signal for each fiber coordinate
% We are not sure about which coordinate is the xyz
% We are not sure how to get the S0 value out of the b=0 (non-diffusion
% weighted) image

%% Please explain

fe_bvecs2             = bvecs'; %feGet(fe,'bvecs');
fe_bvals2             = bvals'./1000; %feGet(fe,'bvals');

%% Comupute prediction diffusion signals
for ii = 1:length(small_fg_img.fibers);
    
    % These coordinates are within the first fiber, and thus the connectome
    % This is an N x 3 matrix of integers, needed by feGet.
    % This way of getting the fiber coordinates matches the coordinates in the
    % ROI.  But it is obviously ridiculous and if we have to have a difference
    % in the coordinate frames, then there has to be a roi2fg() coordinate
    % transform.  We can't write code like this everywhere.
    fCoords = ((uint32(ceil(small_fg_img.fibers{ii})))+1)';
    
    % take S0 image
    S0 = dwiGet(dwi,'b0 image',fCoords); % S0_All = fe.diffusion_S0_img;    
   
    % Compute predicted diffusion direction in each voxel   
    for jj = 1:length(fCoords)
        % For each fiber coordinate create a predicted signal.  Here is the
        % fiber tensor at that coordinate
        voxTensor = Q{ii}(jj,:);
        
        % Add the new signal to the current signals at that coordinate        
        curSig = pData(fCoords(jj,1),fCoords(jj,2),fCoords(jj,3),1:length(fe_bvecs2));
        newSig = wgts(ii)*feComputeSignal(S0(jj), fe_bvecs2, fe_bvals2, voxTensor);
        
        % Store back in
        pData(fCoords(jj,1),fCoords(jj,2),fCoords(jj,3),1:length(fe_bvecs2)) = ...
            squeeze(curSig) + newSig;
        
        % BW to make the visualization work ... also, look over dwiSet/Get/Create
        %
        tmp = feComputeSignal(S0(jj), fe_bvecs2, fe_bvals2, voxTensor);
        dwi2 = dwiSet(dwi,'sig',tmp);
        mrvNewGraphWin; dwiPlot(dwi2,'adc',tmp(11:end),reshape(voxTensor,3,3)) 
        dwiPlot(dwi2,'dsig image xy',tmp)
        
        %           pData(fCoords(jj,1),fCoords(jj,2),fCoords(jj,3),1:length(fe_bvecs2)) =...
        %             wgts(ii)*feComputeSignal(S0(jj), fe_bvecs2, fe_bvals2, voxTensor);
        
        % For the moment we are keeping everything, but we may not have to
        % do this in the long run.
        % OK, keep PSig based on fiber and voxel
        % The empty slots have no fibers.
        % The other slots of 106 diffusion signals predicted
        PSig_voxel{ii,jj} = newSig;

    end
    clear S0
end

%% Franco did get the pSig like this
% Predicted measured signal from both fiber and iso
    %
    % pSig = feGet(fe,'pSig full')
    % pSig = feGet(fe,'pSig full',coords)
    % pSig = feGet(fe,'pSig full',voxelIndices)
    val = [feGet(fe,'M fiber'), feGet(fe,'M iso')] * feGet(fe,'full weights');
   
    %%
    if ~isempty(varargin)
      % voxelIndices     = feGet(fe,'voxelsindices',varargin);
      % voxelRowsToKeep  = feGet(fe,'voxel rows',voxelIndices);
      % val           = val(voxelRowsToKeep,:);
      val = val(feGet(fe,'voxel rows',feGet(fe,'voxelsindices',varargin)));
    end

%% Put pSig in nifti structure
pNifti = niftiSet(pNifti,'data',pData);

% save the nifti file 
if save_flag,
    niftiWrite(pNifti, pNifti.fname)
end

% Let's figure out how to see what we created in the pNifti data set.

return

%% EXTRA REMOVE SOME DAY

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

%% ---------
% Fiber density statistics.
% Computes the fiber density (how many fibers in each voxel)
% We compute the following values:
% (1) The number of fibers in each voxel
% (2) The number of unique fibers with non-zero weights
% (3) The sum of the weights in each voxel
% (4) The mean of the weights in each voxel
% (5) The variance of the weigths in each voxel

FiberStats = feGet(fe,'fiberdensity');






