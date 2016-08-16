function Y = projInfty( X )
%PROJINFTY Summary of this function goes here
%   Detailed explanation goes here
Y = zeros(size(X));
AX = max(sqrt(abs(X(:,:,:,1)).^2 + abs(X(:,:,:,2)).^2 + abs(X(:,:,:,3)).^2),1);
Y(:,:,:,1) = X(:,:,:,1)./AX;
Y(:,:,:,2) = X(:,:,:,2)./AX;
Y(:,:,:,3) = X(:,:,:,3)./AX;

end

