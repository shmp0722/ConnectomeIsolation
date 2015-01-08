function pNifti =  predictDWInifti(fe, dwi)
%
% pNifti =  predictDWInifti(fe, dwi)
%
% This function gives you prediction diffusion nifti file based on LiFE
% weights
% 
% INPUT
% fe  ; After running feConnectomeInit and fitting a model
% dwi ; Structure made by dwiCreate
% 
% Example
% fe;
%   dwiFile       = fullpath to 'dwi1st_b2000_150dirs.nii.gz';
%   fgFileName    = fullpath to 'fgFileName.mat or .pdb'; 
%   feFileName    = 'life_build_model_demo_CSD_PROB';
%   dwiFileRepeat = fullpath to  'dwi2nd_b2000_150dirs.nii.gz';
%   t1File = full path to 't1.nii.gz';
%   % Run LiFE
%   fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],dwiFileRepeat,t1File);
%
%   % fit the model 
%   fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));
% 
% dwi;
%   bvecs =   dlmread('dwi1st_b2000_150dirs.bvecs');
%   bvals =   dlmread('dwi1st_b2000_150dirs.bvals');
%   dwi   = dwiCreate('nifti',dwiFile,'bvecs',bvecs','bvals',bvals');
%
% See; feConnectomeBuildModel, feGet ...
% SO wrote 2015

%% Here is the predicted signal from the life fit for the whole connectome

pSig = feGet(fe,'pSig full');
coords  = feGet(fe,'roi coords');
nVoxels = feGet(fe,'nVoxels');
nBvecs  = feGet(fe,'nbvecs');
nB0     = length(find(dwi.bvals==0));

% Take the pSig values and put them in a data volume like the nifti data
% volume
roiSig = zeros(size(dwi.nifti.data));
for cc=1:nVoxels
    roiSig(coords(cc,1),coords(cc,2),coords(cc,3),(nB0+1):end) =...
        pSig((cc-1)*nBvecs + (1:nBvecs));
end
%% Make an empty nifti file the same size as the original
% check 
% if ~struct(dwi) && ischar(dwi);
%     dwi = niftiRead(dwi);
% end
    
% Put pSig back to original nifti
pData = dwi.nifti.data;
for cc=1:nVoxels
    pData(coords(cc,1),coords(cc,2),coords(cc,3),11:end) = pSig((cc-1)*nBvecs + (1:nBvecs));
end

%% duplicate original nifti structure
pNifti      = dwi.nifti;
pNifti.data = pData;

% strip extension
[p,f] = fileparts(pNifti.fname);
[~,f] = fileparts(f);

% give correct extension to the file
newFname = fullfile(p,[f,'_Psig.nii.gz']);
pNifti.fname = newFname;


