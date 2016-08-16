function [D,Dt] = defDDt
        % defines finite difference operator D
        % and its transpose operator
        
        D = @(U) grad(U);
        Dt = @(V) div(V);
        function DU = grad(U)
            % Forward finite difference operator 
            %(with circular boundary conditions)
            DU(:,:,:,1) = [diff(U,1,2), U(:,1,:) - U(:,end,:)]; % x corresponds to columns (2nd element)
            DU(:,:,:,2) = [diff(U,1,1); U(1,:,:) - U(end,:,:)]; % y corresponds to rows (1st element)
            DU(:,:,:,3) = cat(3, diff(U,1,3), U(:,:,1) - U(:,:,end)); % z corresponds to 3rd element
        end
        
        function DtXYZ = div(V)
            % Divergence operator (transpose of gradient)
            X = V(:,:,:,1);
            Y = V(:,:,:,2);
            Z = V(:,:,:,3);
            DtXYZ = [X(:,end,:) - X(:,1,:), -diff(X,1,2)];
            DtXYZ = DtXYZ + [Y(end,:,:) - Y(1,:,:); -diff(Y,1,1)];
            DtXYZ = DtXYZ + cat(3, Z(:,:,end) - Z(:,:,1), -diff(Z,1,3));
        end
end

