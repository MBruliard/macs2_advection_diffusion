function [ R ] = rigidityMatrix( PDE, Nx, Ny, dx, dy )
% rigidityMatrix: return the rigidity matrix R.
% [ R ] = rigidityMatrix( PDE, Nx, Ny, dx, dy );
%
% Inputs:
%  PDE: PDE object containing information about the system
%  Nx:  number of discretization columns on the x axis
%  Ny:  number of discretization rows on the y axis
%  dx:  step size in the Ox direction
%  dy:  step size in the Oy direction
%
% Outputs:
%  R: matrix of size Ny x Nx by Ny x Nx, where Nx is the number of columns
%     and Ny the number of rows of the mesh


% R = sparse(Ny*Nx, Ny*Nx);

% A/B -> 3*(Nx-2) + 4
% R   -> (3*(Ny-2) + 4) * (3*(Nx-2) + 4)
nb_nonzero = (3*(Ny-2) + 4) * (3*(Nx-2) + 4);
R = spalloc(Ny*Nx, Ny*Nx, nb_nonzero);

for ny = 1 : Ny

    if ny > 1
        % We define the matrix B at position (..., dy*ny)
        for nx = 1 : Nx
            k = nx + (ny-1) * Nx;
            x_pos = PDE.x_min + dx*(nx+1); % Current x position
            y_pos = PDE.y_min + dy*(ny+1); % Current y position
            if nx > 1
                delta =   PDE.K.eval_coef(1,2,x_pos-dx,y_pos)*dx*dy + ...
                        + PDE.K.eval_coef(2,1,x_pos,y_pos-dy)*dx*dy;
                R(k,k-Nx-1) = delta;
            end
            gamma =   PDE.K.eval_coef(1,2,x_pos-dx,y_pos)*dx*dy + ...
                    - PDE.K.eval_coef(1,2,x_pos+dx,y_pos)*dx*dy + ...
                    + 4*PDE.K.eval_coef(2,2,x_pos,y_pos-dy)*dx^2;
            R(k,k-Nx) = gamma;
            if nx < Nx
                delta =   PDE.K.eval_coef(1,2,x_pos+dx,y_pos)*dx*dy + ...
                        + PDE.K.eval_coef(2,1,x_pos,y_pos-dy)*dx*dy;
                R(k,k-Nx+1) = -delta;
            end
        end
    end

    % We define the matrix A at position (..., dy*ny)
    for nx = 1 : Nx
        k = nx + (ny-1) * Nx;
        x_pos = PDE.x_min + dx*(nx+1); % Current x position
        y_pos = PDE.y_min + dy*(ny+1); % Current y position
        if nx > 1
            beta =   PDE.K.eval_coef(2,1,x_pos,y_pos-dy)*dx*dy + ...
                   - PDE.K.eval_coef(2,1,x_pos,y_pos+dy)*dx*dy + ...
                   + 4*PDE.K.eval_coef(1,1,x_pos-dx,y_pos)*dy^2;
            R(k,k-1) = beta;
        end
        alpha = - 4*PDE.K.eval_coef(1,1,x_pos-dx,y_pos)*dy^2 + ...
                - 4*PDE.K.eval_coef(1,1,x_pos+dx,y_pos)*dy^2 + ...
                - 4*PDE.K.eval_coef(2,2,x_pos,y_pos-dy)*dx^2 + ...
                - 4*PDE.K.eval_coef(2,2,x_pos,y_pos+dy)*dx^2;
        R(k,k) = alpha;
        if nx < Nx
            beta =   PDE.K.eval_coef(2,1,x_pos,y_pos+dy)*dx*dy + ...
                   - PDE.K.eval_coef(2,1,x_pos,y_pos-dy)*dx*dy + ...
                   + 4*PDE.K.eval_coef(1,1,x_pos+dx,y_pos)*dy^2;
            R(k,k+1) = beta;
        end
    end

    if ny < Ny
        % We define the matrix B^t at position (..., dy*ny)
        for nx = 1 : Nx
            k = nx + (ny-1) * Nx;
            x_pos = PDE.x_min + dx*(nx+1); % Current x position
            y_pos = PDE.y_min + dy*(ny+1); % Current y position
            if nx > 1
                delta =   PDE.K.eval_coef(1,2,x_pos-dx,y_pos)*dx*dy + ...
                        + PDE.K.eval_coef(2,1,x_pos,y_pos+dy)*dx*dy;
                R(k,k+Nx-1) = -delta;
            end
            gamma =   PDE.K.eval_coef(1,2,x_pos+dx,y_pos)*dx*dy + ...
                    - PDE.K.eval_coef(1,2,x_pos-dx,y_pos)*dx*dy + ...
                    + 4*PDE.K.eval_coef(2,2,x_pos,y_pos+dy)*dx^2;
            R(k,k+Nx) = gamma;
            if nx < Nx
                delta =   PDE.K.eval_coef(1,2,x_pos+dx,y_pos)*dx*dy + ...
                        + PDE.K.eval_coef(2,1,x_pos,y_pos+dy)*dx*dy;
                R(k,k+Nx+1) = delta;
            end
        end
    end

end


end