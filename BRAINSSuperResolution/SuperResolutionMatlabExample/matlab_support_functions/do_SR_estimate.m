function [X_hr] = do_SR_estimate(X_lr,edgemap)
%%
% define the downsampling operator (A).
% Here A is defined by downsampling in fourier domain. However, it is not
% necessarily the best downsampling operator, and any downsampling operator
% can be used.
lowres = size(X_lr);
highres = size(edgemap);
k = get_kspace_inds_pos( highres ); %k=fourier indices
ind_samples = get_lowpass_inds_pos(k,lowres);

%--HJ [A,At] = defAAt_fourier(ind_samples, highres); %Define function handles for fourier projection operators
%% Rescale data -- X_hr will have same scaling as X_lr_pad
FX_lr_pad = complex(zeros(highres,'single'),0);
fftTemp = single(fftn(X_lr));
FX_lr_pad(ind_samples) = single(fftTemp);
X_lr_pad = single((prod(highres)/prod(lowres))*real(ifftn(FX_lr_pad)));

b = single(A_fhp(X_lr_pad,ind_samples, highres));


%% Run Weighted L2 algorithm - ADMM version
lambda = single(1e-3); %regularization parameter
Niter = single(100);  %number of iterations
gam = single(1);   %ADMM parameter, in range [0.1,10]
tol = single(1e-8);  %convergence tolerance
%--HJ tic
[X_hr, cost, resvec ] = OpWeightedL2(b,edgemap,lambda,ind_samples,highres,Niter,tol,gam);
%X_hr = NormalizeDataComponent(abs(X_hr));
%plot(cost)
%plot(resvec)
end
