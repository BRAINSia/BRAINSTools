function Y = shrink(X,tau)
%Returns thresholded magnitudes of X.
Z = abs(X);
Z(Z<eps) = 1;
Y = max(abs(X) - tau, 0)./Z;
end


