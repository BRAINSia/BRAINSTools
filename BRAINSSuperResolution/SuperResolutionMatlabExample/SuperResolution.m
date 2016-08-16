function [ X_sr ] = SuperResolution( edgemap, dwi_b0_lr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% edgemap_size = size(edgemap);
% edgemap_2d = edgemap(:,:,round(edgemap_size(3)/2)); % make edgemask 2D %%%%%%%%%%%%
% figure(1); imagesc(edgemap_2d,[0 0.5]); colorbar; title('spatial weights image');

[X_lr, oldMin, scaleFactor] = NormalizeDataComponent(single(dwi_b0_lr));

% X_lr_size = size(X_lr);
% X_lr_2d = X_lr(:,:,round(X_lr_size(3)/2)); % make edgemask 2D %%%%%%%%%%%%
% figure(2); imagesc(X_lr_2d,[0 1]); colorbar; title('input low resolution image');

%% Run Super-Resolution reconstruction
X_sr = do_SR_estimate(X_lr,edgemap);
X_sr = single(X_sr);
%X_sr = single((X_sr-0.0)*(1.0/scaleFactor)+oldMin);

end

