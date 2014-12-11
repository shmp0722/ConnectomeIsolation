
function [pSig_nifti] = CreatePredictedDiffusionNifti_Voxel(fg, dwiFile,dwiFileRepeat)
% This function returns nifti file of predicted difusion signal
% 
%  [pSig_nifti] = CreatePredictedDiffusionNifti_Voxel(fg, dwi)
% 
% Imput
%  fg = fgRead(fg), fg
%
% SO @ Vista lab, 2014

%% Build the file names for the diffusion data, the anatomical MRI.

% Here are the original dMRI data
dwiFile       = fullfile(lifeDemoDataPath('diffusion'),'life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz');
dwiFileRepeat = fullfile(lifeDemoDataPath('diffusion'),'life_demo_scan2_subject1_b2000_150dirs_stanford.nii.gz');
t1File        = fullfile(lifeDemoDataPath('anatomy'),  'life_demo_anatomy_t1w_stanford.nii.gz');

%% retreive bvals, bvecs from dwi structure
if ischar(dwi),
    dwi = niftiRead(dwi);
elseif isstruct(dwi);
    dwi =dwi;
end
bvecs =   dwi.bvecs;
bvals =   dwi.bvals;
%% Take the first 10 fibers from the fiber group

if ischar(fg),
    fg = niftiRead(fg);
elseif isstruct(fg);
    fg =fg;
end
%% Build a model of the connectome.
%
% Running with the issue2 branch on LiFE

feFileName    = 'tmp_fe';
fe = feConnectomeInit(dwiFile, fg, feFileName,fullfile(fileparts(fgFileName)),dwiFileRepeat,t1File);
%fe = feConnectomeInit(dwiFile, fgFileName, feFileName,fullfile(fileparts(fgFileName)),dwiFileRepeat,t1File);

%% Fit the model with global weights.

fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));
%% the weights from the LiFE solution

% how can we get weights without runnning whole LiFE sequence?
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

% nFibers = length(fe.fg.fibers);
 nFibers =1;

% These coordinates are within the first fiber, and thus the connectome
% This is an N x 3 matrix of integers, needed by feGet.
% This way of getting the fiber coordinates matches the coordinates in the
% ROI.  But it is obviously ridiculous and if we have to have a difference
% in the coordinate frames, then there has to be a roi2fg() coordinate
% transform.  We can't write code like this everywhere.
fCoords = ((uint32(ceil(fe.fg.fibers{nFibers})))+1)';
% cCoords = uint32(feGet(fe,'roi coords'));

% This takes an 
foundVoxels = feGet(fe,'find voxels',fCoords);

% take S0 image 
% S0 = feGet(fe,'b0 signal image'); % All
S0 = dwiGet(dwi,'b0 image',fCoords);
%  S0_All = fe.diffusion_S0_img; 
%  S0_voxel = dwiGet(dwi,'b0 image',fCoords(1,:));


%         % Return an S0 value for each voxel in coords
%         % Coords are in the rows (i.e., nCoords x 3)
%         if ~isempty(varargin), coords = varargin{1};
%         else error('coords required');
%         end
%         coords = coordCheck(coords);
%         
%         % Indices of the coords in the 3D volume.
%         indx = sub2ind(dwi.nifti.dim(1:3),coords(:,1),coords(:,2),coords(:,3));
%         b0 = dwiGet(dwi,'b0 image nums');
%         val = zeros(size(coords,1),length(b0));
%         for ii = 1:length(b0)
%             tmp = squeeze(dwi.nifti.data(:,:,:,b0(ii)));
%             val(:,ii) = tmp(indx);
%         end

%     Q{nFinbers}(nNodes,:);

fe_bvecs             = feGet(fe,'bvecs');
fe_bvals             = feGet(fe,'bvals');

for jj = 1:length(foundVoxels)
    % Once we get the S0 values for this particular voxel, we can compute
   
    voxTesorQ = Q{nFibers}(jj,:);
    % voxTesorQ = feGet(fe,'voxel tensors',foundVoxels(jj));
    % TensorVoxel = dwiGet(dwi,'tensor image', fCoords(jj,:)); ??

    voxPSig{nFibers,jj} = feComputeSignal(S0(jj), fe_bvecs, fe_bvals, voxTesorQ);
end

% Size; S0 [1,1], bvecs[96, 3], bvals[96, 1] , Q[2,9]
% fe_bvals and S0 are still unclear  
nifti.data(73,29,33,:)

% for ii = 1:length(b0)
%     tmp = squeeze(dwi.nifti.data(:,:,:,b0(ii)));
%     val(:,ii) = tmp(indx);
% end



%%
for ii=1:length(foundVoxels)
    for jj=1:length(fe_bvecs)
        pData(oneFiber(1,1),oneFiber(2,1),oneFiber(3,1),jj) = ...
            wgts(ii)*feComputeSignal(S0, bvecs', bvals(:), Q{1}); 
    end
end

pNifti = niftiSet(pNifti,'data',pData);


%%
for ii=1:length(foundVoxels)
    for jj=1:length(bvecs)
        pData(oneFiber(1,1),oneFiber(2,1),oneFiber(3,1),jj) = ...
            wgts(ii)*feComputeSignal(S0, bvecs', bvals(:), Q{1}); 
    end
end

pNifti = niftiSet(pNifti,'data',pData);

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


%---------
% Given an VOI finds the indices of the matching voxels inside the big
% volume, which ordinarily represents the full connectome.
%
foundVoxels = feGet(fe,'find voxels',coords);
%
% coords is a Nx3 set of coordinates in image space
% foundVoxels is a vector of 1's and 0's, there is a one for each
% location the the connectome coordinates for which there is a match in
% the coords


% We are not sure about which coordinate is the xyz
% We are not sure how to get the S0 value out of the b=0 (non-diffusion
% weighted) image
oneFiber = floor(fe.fg.fibers{1});

VOI_coords = feGet(fe,'roi coords');


% We want the S0 from the raw data, and then we want the S0 values for each
% voxel in the fiber
S0 = feGet(fe,'b0 signal image');
% feGet(fe,'voxels indices',fe.fg.fibers)

% val = feGet(fe,'b0 signal image',int32(oneFiber));
  
% Once we get the S0 values for this particular voxel, we can compute
voxPSig = feComputeSignal(S0(1), bvecs', bvals(:), Q{1});

for ii=1:length(oneFiber)
    for jj=1:length(bvec)
        pData(oneFiber(1,1),oneFiber(2,1),oneFiber(3,1),jj) = ...
            wgts(ii)*exp(-b*(bvec(jj)'*Q*bvec(jj))); 
    end
end


pNifti = niftiSet(pNifti,'data',pData);
%%
% This returns a matrix that is size nBvecs x num_unique_fibers 
  voxelPSignal      = feComputeVoxelSignal(fe,1);
  voxTensors        = feGet(fe,'voxeltensors', 1);  % Get the tensors for each node in each fiber 

% Index for the voxel
%     vv = feGet(fe,'voxelsindices');
    
    nVoxels = feGet(fe,'n voxels'); % Number of voxels in the connectome/VOI
    
    
    val = fe.diffusion_S0_img(feGet(fe,'voxelsindices',nVoxels), :);
    
    vv = 1:nVoxels;
    usedVoxels = feGet(fe,'usedVoxels');
    voxIndex   = usedVoxels(vv);
    nNodes     = feGet(fe,'nNodes');
    function CreatePredictedDiffusionNifti_Voxel
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

% clear bvecs, clear bvals;
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
% the weights from the LiFE solution
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
% each of its neighbors
% The diagonal parameters are for a perfect stick, [1 0 0]
% Q = feComputeCanonicalDiffusion(fe.fg.fibers, [1 0 0]); 
Q = feGet(fe,'tensors');

%% Just thinking about in VOI
% We don't need the Q value in each fiber and nodes

% Return the VOI comprised by the full connectome.
VOI_coords = feGet(fe,'roi coords'); % in img space
%
% S0 = feGet(fe,'b0 signal image',VOI_coords(1,:)');%,coords');
S0 = feGet(fe,'b0 signal image');%,VOI_coords(1,:)');%,coords');

% get Diffusion direction in each voxels
Q2 = feComputeCanonicalDiffusion_voxel(VOI_coords', [1 0 0]); % Q =voxTensors;

for ii = 1:length(VOI_coords)
    voxDSig{ii} = feComputeSignal(S0, bvecs', bvals(:), Q);
end

% for ii=1:length(VOI_coords)
    for jj=1:length(bvecs)
        pData(VOI_coords(1,1),VOI_coords(2,1),VOI_coords(3,1),jj) = ...
            feGet(fe,'pSig f vox', VOI_coords(jj,:)');

    end
% end


%% Add diffusion signal for each fiber coordinate

% We are not sure about which coordinate is the xyz
% We are not sure how to get the S0 value out of the b=0 (non-diffusion
% weighted) image
oneFiber = VOI_coords';

% We want the S0 from the raw data, and then we want the S0 values for each
% voxel in the fiber
% S0 = feGet(fe,'b0 signal image');
% feGet(fe,'voxels indices',fe.fg.fibers)

% val = feGet(fe,'b0 signal image',int32(oneFiber));
  
% Once we get the S0 values for this particular voxel, we can compute
for ii = 1:length(oneFiber)
    voxDSig{ii} = feComputeSignal(S0(ii), bvecs, bvals, Q2);
end


for ii=1:length(oneFiber)
    for jj=1:length(bvecs)
        pData(oneFiber(1,1),oneFiber(2,1),oneFiber(3,1),jj) = ...
            wgts(ii)*exp(-b*(bvec(jj)'*Q*bvec(jj))); 
    end
end

pNifti = niftiSet(pNifti,'data',pData);


%% Add diffusion signal for each fiber coordinate
% Original; 

% We are not sure about which coordinate is the xyz
% We are not sure how to get the S0 value out of the b=0 (non-diffusion
% weighted) image
oneFiber = floor(fe.fg.fibers{1});

% We want the S0 from the raw data, and then we want the S0 values for each
% voxel in the fiber
% S0 = feGet(fe,'b0 signal image');
feGet(fe,'voxels indices',fe.fg.fibers)

val = feGet(fe,'b0 signal image',int32(oneFiber));
  
% Once we get the S0 values for this particular voxel, we can compute
voxDSig = feComputeSignal(S0, bvecs', bvals(:), Q{1});

for ii=1:length(oneFiber)
    for jj=1:length(bvec)
        pData(oneFiber(1,1),oneFiber(2,1),oneFiber(3,1),jj) = ...
            wgts(ii)*exp(-b*(bvec(jj)'*Q*bvec(jj))); 
    end
end


pNifti = niftiSet(pNifti,'data',pData);

%%
% This returns a matrix that is size nBvecs x num_unique_fibers 
  voxelPSignal      = feComputeVoxelSignal(fe,1);
  voxTensors        = feGet(fe,'voxeltensors',1);  % Get the tensors for each node in each fiber 

% Index for the voxel
%     vv = feGet(fe,'voxelsindices');
    
    nVoxels = feGet(fe,'n voxels'); % Number of voxels in the connectome/VOI    
    
    val = fe.diffusion_S0_img(feGet(fe,'voxelsindices',nVoxels), :);
    
    vv = 1:nVoxels;
    usedVoxels = feGet(fe,'usedVoxels');
    voxIndex   = usedVoxels(vv);
    nNodes     = feGet(fe,'nNodes');
    
     % Get the tensors for each node in each fiber going through this voxel:
    val = zeros(nNodes(voxIndex), 9); % Allocate space for all the tensors (9 is for the 3 x 3 tensor components)
    for ii = 1:nNodes(voxIndex)           % Get the tensors
      val(ii,:) = fe.fibers.tensors{fe.voxel2FNpair{voxIndex}(ii,1)} ...
        (fe.voxel2FNpair{voxIndex}(ii,2),:);
    end

%% 
feFileName    = 'life_build_model_demo_CSD_PROB_small';

%% dwi Get
% Keep original index
indx = sub2ind(dwi.nifti.dim(1:3),VOI_coords(:,1),VOI_coords(:,2),VOI_coords(:,3));
%% Get diffusion data from a fiber
dSig    = dwiGet(dwi,'diffusion data image',VOI_coords(1,:));
nifti.data(indx(1))

%% Compute the predicted signal by each tensors of each node in this voxel.
usedVoxels  =1;

S0                = feGet(fe,'b0signalimage');
bvecs             = feGet(fe,'bvecs');                      % bvecs
bvals             = feGet(fe,'bvals');                      % bvals
tot_fibers_num    = feGet(fe,'tot f num');  % number of total fibers in the voxel
unique_fibers_num = feGet(fe,'unique f num',    usedVoxels);  % number of unique fibers in the voxel
tot_fiber_index   = cell2mat(feGet(fe,'totf',   usedVoxels)); % indexes to the total fibers in the voxels
unique_fiber_index= cell2mat(feGet(fe,'uniquef',usedVoxels)); % indexes to the unique fibers in the voxels
voxTensors        = feGet(fe,'voxeltensors',usedVoxels);

Q = feComputeCanonicalDiffusion(small_fg.fibers, [1 0 0]); % Q =voxTensors;

voxDSig = feComputeSignal(S0, bvecs, bvals, Q);

%%
% Indexes of actually used voxels, this can be less than the number of
% voxels int he ROI in case some voxels have no fibers in them
usedVoxels   = feGet(fe,'usedVoxels');
nVoxels      = length(usedVoxels);
nBvecs       = feGet(fe,'nBvecs');
vox_sparse_pSig = cell(nVoxels,1);

% % Return the model (M matrix), or a subset of it.
% Mfiber = feGet(fe,'model',coords');
% 
% % Isotropic portion of M matrix. Returns the matrix for the full model or
% % for a subset of voxels.
% 
% Miso = feGet(fe,'M iso',coords');
% 
% % Return the number of the fibers gooes through the voxels.
% idxFibers = feGet(fe,'uniquef',coords');    % for some the voxels,
% % Return the total number of fibers for all the voxels or in a set of
% % voxels
% nFibers = feGet(fe,'totfnum');
%%
% % predicted = feGet(fe,'pSig f vox');
% pSig = feGet(fe,'pSig full', VOI_coords(1,:));
% measured = feGet(fe,'dSig full by Voxel', VOI_coords(1,:)');
% 
% predicted2 = feGet(fe,'pSig f vox', VOI_coords(2,:)');
% measured2 = feGet(fe,'dSig full by Voxel', VOI_coords(2,:)');
% 
% 
% size(pSig)

     % Get the tensors for each node in each fiber going through this voxel:
    val = zeros(nNodes(voxIndex), 9); % Allocate space for all the tensors (9 is for the 3 x 3 tensor components)
    for ii = 1:nNodes(voxIndex)           % Get the tensors
      val(ii,:) = fe.fibers.tensors{fe.voxel2FNpair{voxIndex}(ii,1)} ...
        (fe.voxel2FNpair{voxIndex}(ii,2),:);
    end

%% 
feFileName    = 'life_build_model_demo_CSD_PROB_small';

%% dwi Get
% Take all voxels 
coords = horzcat(small_fg.fibers{:});
% transform in img space 
coords = unique(floor(mrAnatXformCoords(dwi.nifti.qto_ijk,coords)),'rows');
% Keep original index
indx = sub2ind(dwi.nifti.dim(1:3),coords(:,1),coords(:,2),coords(:,3));
%% Get diffusion data from a fiber
dSig    = dwiGet(dwi,'diffusion data image',coords);

%% Compute the predicted signal by each tensors of each node in this voxel.

S0                = feGet(fe,'b0signalimage',usedVoxels);
bvecs             = feGet(fe,'bvecs');                      % bvecs
bvals             = feGet(fe,'bvals');                      % bvals
tot_fibers_num    = feGet(fe,'tot f num');  % number of total fibers in the voxel
unique_fibers_num = feGet(fe,'unique f num',    usedVoxels);  % number of unique fibers in the voxel
tot_fiber_index   = cell2mat(feGet(fe,'totf',   usedVoxels)); % indexes to the total fibers in the voxels
unique_fiber_index= cell2mat(feGet(fe,'uniquef',usedVoxels)); % indexes to the unique fibers in the voxels
voxTensors        = feGet(fe,'voxeltensors',usedVoxels);

Q = feComputeCanonicalDiffusion(small_fg.fibers, [1 0 0]); % Q =voxTensors;

voxPSig = feComputeSignal(S0, bvecs, bvals, Q);

%%
% Indexes of actually used voxels, this can be less than the number of
% voxels int he ROI in case some voxels have no fibers in them
usedVoxels   = feGet(fe,'usedVoxels');
nVoxels      = length(usedVoxels);
nBvecs       = feGet(fe,'nBvecs');
vox_sparse_pSig = cell(nVoxels,1);


