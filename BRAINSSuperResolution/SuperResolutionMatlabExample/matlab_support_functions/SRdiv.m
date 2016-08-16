function DtXYZ = SRdiv(V)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

            % Divergence operator (transpose of gradient)
            X_contiguous= V(:,:,:,2);
            Y_lines =     V(:,:,:,1);
            Z_slices =    V(:,:,:,3);
            DtXYZ_contiguous = [X_contiguous(:,end,:) - X_contiguous(:,1,:), -diff(X_contiguous,1,2)];
            DtXYZ_1 = [Y_lines(end,:,:) - Y_lines(1,:,:); -diff(Y_lines,1,1)];
            DtXYZ_3 = cat(3, Z_slices(:,:,end) - Z_slices(:,:,1), -diff(Z_slices,1,3));
            DtXYZ = DtXYZ_contiguous + DtXYZ_1 + DtXYZ_3;
end
