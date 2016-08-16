function DU = SRgrad(U)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

            % Forward finite difference operator 
            %(with circular boundary conditions)
            DU_contiguous =   [diff(U,1,2), U(:,1,:) - U(:,end,:)]; % x corresponds to columns (2nd element)
            DU_lines =        [diff(U,1,1); U(1,:,:) - U(end,:,:)]; % y corresponds to rows (1st element)
            DU_slices = cat(3, diff(U,1,3), U(:,:,1) - U(:,:,end)); % z corresponds to 3rd element
            
            DU = cat(4, DU_lines, DU_contiguous, DU_slices);
end

