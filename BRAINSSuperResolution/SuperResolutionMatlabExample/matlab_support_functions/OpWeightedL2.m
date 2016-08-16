function [X, cost, resvec] = OpWeightedL2(b_FC,edgemask,lambda,ind_samples,res,iter,tol,gam)
% OPWEIGHTEDL2: Solves weighted L2 regularized inverse problems.
% Minimizes the cost function
% X* = argmin_X ||A(X)-b_FC||_2^2 + lambda ||W |D(X)| ||_2^2
% where     X* = recovered image
%           A  = linear measurement operator
%           b_FC  = (noisy) measurements
%           W  = diagonal weight matrix built from the edge mask
%           |D(X)| = gradient magnitude at each pixel
%
% Inputs:  A = function handle representing the forward
%               model/measurement operator
%          At = function handle representing the backwards model/
%               the transpose of the measurment operator.
%               (e.g. if A is a downsampling, At is a upsampling)
%          b_FC =  a vector of measurements; should match the
%               dimensions of A(X)
%          lambda = regularization parameter that balances data fidelity
%               and smoothness. set lambda high for more smoothing.
%          siz = output image size, e.g. siz = [512,512]
%          Niter = is the number of iterations; should be ~100-500
%
% Output:  X = high-resolution output image
%          cost = array of cost function value vs. iteration
%Define AtA fourier mask
p_image = zeros(res,'single'); 
p_image(1,1,1) = single(1);

AtA =single( At_fhp(single(A_fhp(p_image,ind_samples,res)),ind_samples,res));
AtAhat = single(fftn(AtA));

divImg = SRdiv(SRgrad(p_image));

grad_p_img = SRgrad(p_image);

DtDhat = single(fftn(SRdiv(SRgrad(p_image))));

mu = single(edgemask);
Atb = single(At_fhp(b_FC,ind_samples,res));
X = Atb;

DX = SRgrad(X);

L = single(zeros(size(DX)));
resvec = single(zeros(1,iter));
cost = single(zeros(1,iter));
TwoAtb = single(2)*Atb;
twoMuPlusGam = single((2*mu + gam).^(-1));
template.data = single(twoMuPlusGam);

muinv = repmat(twoMuPlusGam,[1,1,1,3]);
denominator = (2*AtAhat + lambda*gam*DtDhat);

for i = 1:iter
    disp(i)
    % Y subprob
    Z = gam*(DX+L);
    Y = muinv.*Z;

    % X subprob
    numerator = fftn(TwoAtb + lambda*gam*SRdiv(Y-L));
    temp = numel(numerator)*ifftn(numerator);
    X = single(real(ifftn(numerator./denominator)));
    % L update
    DX = SRgrad(X);
    residue = DX-Y;
    resvec(i) = norm(residue(:))/norm(Y(:));
    if (i > 10) && (resvec(i) < tol)
        return;
    end
    L = L + residue;
    %Calculate cost function
    WDX = repmat(sqrt(mu),[1,1,1,3]).*DX;
    diff = A_fhp(X,ind_samples,res)-b_FC;
    cost(i) = norm(diff(:)).^2 + lambda*norm(WDX(:)).^2; %cost function
end
end
