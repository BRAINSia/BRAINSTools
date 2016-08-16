function [ out ] = A_fhp( z, ind_samples, res )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        p = 1/sqrt(prod(res))*fftn(z);
        out = p(ind_samples);
end

