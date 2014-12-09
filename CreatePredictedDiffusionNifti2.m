function CreatePredictedDiffusionNifti2
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

%% Build a model of the connectome.
%
% Running with the issue2 branch on LiFE
%
% If we don't need fiber weights, we can skip this

feFileName    = 'life_build_model_demo_CSD_PROB_small';
fe = feConnectomeInit(dwiFile, small_fg, feFileName,fullfile(fileparts(fgFileName)),dwiFileRepeat,t1File);

%% Fit the model with global weights.

fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));

%% the weights from the LiFE solution

wgts = feGet(fe,'fiber weights');

%% Take the first 10 fibers from the fiber group

% This is intended to run quickly
fgFileName    = fullfile(lifeDemoDataPath('tractography'), ...
    'life_demo_mrtrix_csd_lmax10_probabilistic.mat');

% pick first 10 fibers as nwe fg
fg = fgRead(fgFileName);
small_fg = fgCreate('name', ['small_' fg.name], 'fibers', fg.fibers(1:10));
small_fg.pathwayInfo = fg.pathwayInfo(1:10);
clear fg

% transform coordinates in image space form acpc
small_fg_img = dtiXformFiberCoords(small_fg, inv(dwi.nifti.qto_xyz), 'img');

% gather all img voxels
img_coords  = horzcat(small_fg_img.fibers{:});
img_coords  = (ceil(img_coords))'+1;

% unique_coords = unique(img_coords,'rows');
% Get voxcel indices
[~,voxIndices] = ismember(img_coords, feGet(fe,'roi coords'),'rows'); % This is slow
% [ind,voxIndices] = ismember(unique_coords, feGet(fe,'roi coords'),'rows'); % This is slow

% %%
% % Indexes of actually used voxels, this can be less than the number of
% % voxels int he ROI in case some voxels have no fibers in them
% % voxIndices   = feGet(fe,'voxIndices');
% nVoxels      = length(voxIndices);
% nBvecs       = feGet(fe,'nBvecs');
% vox_sparse_pSig = cell(nVoxels,1);
%
% % Generate the signal predicted in each voxel by the connectome.
% %
% % This for-loop is written in such a way that matlab (ver later than 7.9)
% % will run it in parallel if parallel processing is allowed.
% fprintf('LiFE - Predicting diffusion signal in %i voxel...\n',nVoxels);
%
% % feOpenLocalCluster
% for vv = 1:nVoxels
%   num_unique_fibers = feGet(fe,'unique f num',voxIndices(vv));
%   % This returns a matrix that is size nBvecs x num_unique_fibers
%   voxelPSignal(1:nBvecs,vv) = feComputeVoxelSignal(fe,voxIndices(vv));
% end

%% for example vv = 88, the voxel has two unique fiber
for  vv =1: length(voxIndices)
    % Extract information regarding, voxels, signal and fibers.
    S0                = feGet(fe,'b0signalimage',   voxIndices(vv));  % non diffusion-weighted signal
    bvecs             = feGet(fe,'bvecs');                      % bvecs
    bvals             = feGet(fe,'bvals');                      % bvals
    tot_fibers_num    = feGet(fe,'tot f num',       voxIndices(vv));  % number of total fibers in the voxel
    unique_fibers_num = feGet(fe,'unique f num',    voxIndices(vv));  % number of unique fibers in the voxel
    tot_fiber_index   = cell2mat(feGet(fe,'totf',   voxIndices(vv))); % indexes to the total fibers in the voxels
    unique_fiber_index= cell2mat(feGet(fe,'uniquef',voxIndices(vv))); % indexes to the unique fibers in the voxels
    voxTensors        = feGet(fe,'voxeltensors',    voxIndices(vv));  % Get the tensors for each node in each fiber
    
    
    % Compute the predicted signal by each tensors of each node in this voxel.
    voxDSig = feComputeSignal(S0, bvecs, bvals, voxTensors);
    
    % We (FRK/BW) now think we should never really enter this condition.
    % At some point there were redundant representations of the nodes, and they
    % could be eliminated here.  But we don't think this happens any more.
    %
    % The code here found the redundant nodes and removed them by
    % multiplying with the combineM matrix.
    % Use only the prediction from the unique fibers, not from all the fibers.
    if tot_fibers_num ~= unique_fibers_num
        warning('redundant nodes')
        combineM = zeros(tot_fibers_num, unique_fibers_num);
        
        % The matrix combineM is a set of 0s and 1s that will sum together the
        % nodes from a single fiber.
        for ii=1:unique_fibers_num
            combineM(:,ii) = ( tot_fiber_index == unique_fiber_index(ii) );
        end
        
        % The matrix for this voxel starts with each node, and when we multiply
        % by combineM. The resulting matrix represents each fiber
        % (not each node) as a column
        VoxDSig(1:nBvecs,vv) = voxDSig*combineM;
        clear voxDSig
    end
end

%   % Fibers in the connectome determine the directional signal in the voxel
%   % signal, not the mean signal in the voxel. Typically, we demean the
%   % voxel signal we will predict.
%   %
%   %  demeaned_pSig = voxelPSignal - repmat(mean(voxelPSignal, 1),nBvecs,1);
%   %
%   % The mean is predicted by the Miso part of the matrix, not the
%   % Mfiber part.
%   %
%   % Then we reorganize the demeaned signal into a column vector
%   %
%   %  vox_sparse_pSig{vv} = reshape(demeaned_pSig', num_unique_fibers * nBvecs, 1) ;
%   %
%   % Somehow this column vector ends up occupying the right parts of the
%   % Mfiber matrix when we are done.  That miracle happens below.
%
%       % In typical application, we remove the mean from the fiber
%       % predictions.  The returned prediction has zero mean.
%       vox_sparse_pSig{vv}   = reshape((voxelPSignal - repmat(mean(voxelPSignal, 1),nBvecs,1))', ...
%           num_unique_fibers * nBvecs, 1) ;
% end
%
%





%% Build a model of the connectome.
%
% Running with the issue2 branch on LiFE
%
% If we don't need fiber weights, we can skip this

feFileName    = 'life_build_model_demo_CSD_PROB_small';
fe = feConnectomeInit(dwiFile, small_fg, feFileName,fullfile(fileparts(fgFileName)),dwiFileRepeat,t1File);

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
newFname = [p,f,'_Predicted.nii.gz'];
% give name to pNifti
pNifti.fname = newFname;
%% Compute the diffusion tensors for each node in each fiber

% The direction at the node is the average direction between this node and
% each of its neighbors% Isotropic portion of M matrix. Returns the matrix for the full model or

% The diagonal parameters are for a perfect stick, [1 0 0]
Q = feComputeCanonicalDiffusion(small_fg_img, [1 0 0]);  % Q = feGet(fe,''fibers tensors');
% Q = feComputeCanonicalDiffusion_voxel(img_coords, [1 0 0]);
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
    
    %% Compute predicted diffusion direction
    % i don't know how fiber weights = wgts(ii)* work (SO)
    for jj = 1:length(fCoords)
        voxTesor = Q{ii}(jj,:);
        pData(fCoords(jj,1),fCoords(jj,2),fCoords(jj,3),1:length(fe_bvecs2)) =...
            feComputeSignal(S0(jj), fe_bvecs2, fe_bvals2, voxTesor);
    end
    clear S0, clear indx
end

% %% weights
% for ii = 1:length(fe.fg.fibers);
%
%     % These coordinates are within the first fiber, and thus the connectome
%     % This is an N x 3 matrix of integers, needed by feGet.
%     % This way of getting the fiber coordinates matches the coordinates in the
%     % ROI.  But it is obviously ridiculous and if we have to have a difference
%     % in the coordinate frames, then there has to be a roi2fg() coordinate
%     % transform.  We can't write code like this everywhere.
%     fCoords = ((uint32(ceil(fe.fg.fibers{ii})))+1)';
%
%     % take S0 image
%     S0 = dwiGet(dwi,'b0 image',fCoords); % S0_All = fe.diffusion_S0_img;
%
%     %% Compute predicted diffusion direction
%
%     for jj = 1:length(fCoords)
%         voxTesor = Q{ii}(jj,:);
%         pData(fCoords(jj,1),fCoords(jj,2),fCoords(jj,3),1:length(fe_bvecs2)) =...
%             wgts(ii)*feComputeSignal(S0(jj), fe_bvecs2, fe_bvals2, voxTesor);
%     end
%     clear S0, clear indx
% end

%% Put pSig in nifti structure
pNifti = niftiSet(pNifti,'data',pData);

return

% save the nifti file
if save_flag,
    niftiWrite(pNifti, pNifti.fname)
end
