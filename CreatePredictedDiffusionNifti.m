function CreatePredictedDiffusionNifti
%% Use Franco's example data

try, lifeDemoDataPath
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

clear bvecs, clear bvals;
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

%% We need the weights from the LiFE solution


%% Build a model of the connectome.
%
% Running with the issue2 branch on LiFE

feFileName    = 'life_build_model_demo_CSD_PROB_small';
fe = feConnectomeInit(dwiFile, small_fg, feFileName,fullfile(fileparts(fgFileName)),dwiFileRepeat,t1File);
%fe = feConnectomeInit(dwiFile, fgFileName, feFileName,fullfile(fileparts(fgFileName)),dwiFileRepeat,t1File);

%% Fit the model with global weights.
fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));

wgts = feGet(fe,'fiber weights');
coords = 


%% Make an empty nifti file the same size as the original
pNifti = niftiCreate;
pData = zeros(size(nifti.data));

%% Compute the diffusion tensors for each node in each fiber

Q = feComputeCanonicalDiffusion(fe.fg.fibers, [1 0 0]); % Q =voxTensors;


%% Add diffusion signal for each fiber coordinate
oneFiber = floor(fe.fg.fibers{1});

for ii=1:length(oneFiber)
    for jj=1:length(bvec)
        pData(oneFiber(1,1),oneFiber(2,1),oneFiber(3,1),jj) = ...
            wgts(ii)*exp(-b*(bvec(jj)'*Q*bvec(jj))); 
    end
end


pNifti = niftiSet(pNifti,'data',pData);

%% For each 

%% 
feFileName    = 'life_build_model_demo_CSD_PROB_small';

%% dwi Get

% Examples:
%   To get diffusion data from a fiber
%   nifti = niftiRead('raw/DTI__aligned_trilin.nii.gz');
%   bvecs =   dlmread('raw/DTI__aligned_trilin.bvecs');
%   bvals =   dlmread('raw/DTI_aligned_trilin.bvals');
%   dwi   = dwiCreate('nifti',nifti,'bvecs',bvecs,'bvals',bvals);
%   coords = [64 64 30;64 64 31; 64 64 32];
%
%   dws = dwiGet(dwi,'diffusion data acpc',coords);
%
%   ADC = dwiGet(dwi,'adc data image',coords);
%
%   SNR = dwiGet(dwi,'b0 snr',coords);
%
% See also:  dwiCreate, dwiSet, dtiGet, dtiSet, dtiCreate


% Take all voxels 
coords = horzcat(small_fg.fibers{:});
% transform in img space 
coords = unique(floor(mrAnatXformCoords(dwi.nifti.qto_ijk,coords)),'rows');
% Keep original index
indx = sub2ind(dwi.nifti.dim(1:3),coords(:,1),coords(:,2),coords(:,3));
%% Get diffusion data from a fiber
dSig    = dwiGet(dwi,'diffusion data image',coords);

%% Compute the predicted signal by each tensors of each node in this voxel.

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




