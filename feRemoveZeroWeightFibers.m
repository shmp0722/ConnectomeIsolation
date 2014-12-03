function [fg_removed, keep] =  feRemoveZeroWeightFibers(fe)
% 
% [fg_removed, keep] =  feRemoveZeroWeightFibers(fe)
%
% After runnnig LiFE. Exclude zero weight fibers.
% 
%
% SO 2014 wrote

 fg = feGet(fe, 'fibers acpc');
 fweight = feGet(fe, 'fiber weights');            
 fweight(fweight>0) = 1;
 %
 keep = fweight;
 fg_removed = fgExtract(fg, logical(fweight), 'keep');
