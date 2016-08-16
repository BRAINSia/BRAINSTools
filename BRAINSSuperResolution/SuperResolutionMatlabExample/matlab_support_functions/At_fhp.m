function [ out ] = At_fhp(z,ind_samples, res)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        p = complex(zeros(res,'single'),0);
        p(ind_samples) = z;
        out = sqrt(prod(res))*real(ifftn(p));
end

