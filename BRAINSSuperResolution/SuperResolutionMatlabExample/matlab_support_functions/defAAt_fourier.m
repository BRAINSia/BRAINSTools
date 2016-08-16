%Defines Fourier Undersampling operators A, At for Compressed Sensing MRI
%Inputs: ind_samples = indicies of sampling locations, res = resolution
%Outputs: A = measurement operator, At = measurement operator transpose
function [A,At] = defAAt_fourier(ind_samples,res)
%--HJ    A = @(z) A_fhp(z,ind_samples,res);
%--HJ    At = @(z) At_fhp(z,ind_samples,res);
%#codegen
A = 0;
coder.extrinsic('A_fhp')
A = A_fhp(ind_samples,res);

At = 0;
coder.extrinsic('At_fhp')
At = At_fhp(ind_samples,res);

end

function out = A_fhp(z, ind_samples, res)
        p = 1/sqrt(res(1)*res(2)*res(3))*fftn(z);
        out = p(ind_samples);
end

function out = At_fhp(z,ind_samples, res)
        p = zeros(res,'single');
        p(ind_samples) = z;
        out = sqrt(res(1)*res(2)*res(3))*real(ifftn(p));
end
