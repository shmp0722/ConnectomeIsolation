function predictDWI
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

%% Take the first nKept fibers from the fiber group
nKept = 75;

% This is intended to run quickly
fgFileName    = fullfile(lifeDemoDataPath('tractography'), ...
    'life_demo_mrtrix_csd_lmax10_probabilistic.mat');

fg = fgRead(fgFileName);
small_fg = fgCreate('name', ['small_' fg.name], 'fibers', fg.fibers(1:nKept));
small_fg.pathwayInfo = fg.pathwayInfo(1:nKept);
% fgWrite(small_fg, fullfile(lifeDemoDataPath('diffusion'), 'small_fg.mat'));

clear fg

%% Build a model of the connectome.
%
% Running with the issue2 branch on LiFE
%
% If we don't need fiber weights, we can skip this
%
% We are getting warnings here ... try to figure this out with Franco.

feFileName    = 'life_build_model_demo_CSD_PROB_small';
fe = feConnectomeInit(dwiFile, small_fg, feFileName,fullfile(fileparts(fgFileName)),dwiFileRepeat,t1File);

%% Fit the model with global weights.

fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));

%% Here is the predicted signal from the life fit for the whole connectome

Mfiber = feGet(fe,'M fiber'); 
Miso   = feGet(fe,'M iso');
wgts   = feGet(fe,'full weights');
pSig   = [Mfiber,Miso]*wgts; % pSig = feGet(fe,'pSig full');

coords  = feGet(fe,'roi coords');
nVoxels = size(coords,1);
nBvecs  = feGet(fe,'nbvecs');
nB0     = length(find(dwi.bvals==0));

% Take the pSig values and put them in a data volume like the nifti data
% volume
roiSig = zeros(size(nifti.data));
for cc=1:nVoxels
    roiSig(coords(cc,1),coords(cc,2),coords(cc,3),(nB0+1):end) =...
        pSig((cc-1)*nBvecs + (1:nBvecs));
end

% Or, to compare the predicted and observed signal, get the observed signal
% into a vector like pSig
sig  = niftiGet(nifti,'data');
oSig = zeros(size(pSig));
for cc = 1:nVoxels
    oSig((cc-1)*96 + (1:96)) = sig(coords(cc,1),coords(cc,2),coords(cc,3),(nB0+1):end);
end

%% Debug, what is wrong with the coordinates if anything ...

mrvNewGraphWin;
plot(oSig(:),pSig(:),'.')
identityLine
xlabel('Measured'); ylabel('Predicted');
title('All fibers')

%% Figure out how to zero out some of the fibers

newWgts = wgts;
nFibers = size(Mfiber,2);
fRemove = 0.50;   % Fraction of fibers to remove (randomly selected)

% As skip gets smaller, we zero out more fibers
nRemove = round(nFibers*fRemove);
% Find a non-statistics toolbox version of randsample.
% Psychtoolbox has one, I think.
lst = sort(randsample(nFibers,nRemove));   
newWgts(lst) = 0;
pSig2   = [Mfiber,Miso]*newWgts;

mrvNewGraphWin;
plot(oSig(:),pSig2(:),'.')
identityLine
xlabel('Measured'); ylabel('Predicted');
title(sprintf('Removed %i out of %i\n',nRemove,nFibers));

%% SO figure out how to zero out some of the fibers

fiberWgts = feGet(fe,'fiber weights');
isoWgts   = feGet(fe,'iso weights');

% remove fibers by percrntile of weights
% In oreder to see the effect, remove high weighted fibers 
NonZeroFwgts = fiberWgts(fiberWgts>0);
prct = 80;
CutOff = prctile(NonZeroFwgts,prct);
newFwgts = fiberWgts;
newFwgts(newFwgts>CutOff) =0;

% combine new fiber weights with iso weights 
newWgts = [newFwgts;isoWgts];

pSig2   = [Mfiber,Miso]*newWgts;

mrvNewGraphWin;
plot(oSig(:),pSig2(:),'.')
identityLine
xlabel('Measured'); ylabel('Predicted');
% title(sprintf('Removed lowest %i out of %i\n',length(newFwgts(newFwgts==0)),nFibers));
title(sprintf('Removed highest %i out of %i\n',length(newFwgts(newFwgts==0)),nFibers));

%% Make an empty nifti file the same size as the original
% Put pSig back to original nifti
pData = nifti.data;
for cc=1:nVoxels
    pData(coords(cc,1),coords(cc,2),coords(cc,3),11:end) = pSig((cc-1)*96 + (1:96));
end
% duplicate original nifti structure
pNifti      = nifti;
pNifti.data = pData;

% strip .extension
[p,f] = fileparts(pNifti.fname);
[~,f] = fileparts(f);

% give name to pNifti
newFname = [p,f,'_Predicted.nii.gz'];
pNifti.fname = newFname;

%%
return

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






