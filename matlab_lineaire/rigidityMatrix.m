function [ R ] = rigidityMatrix( PDE, Nx, Ny, dx, dy )
% rigidityMatrix: return the rigidity matrix R.
% [ R ] = rigidityMatrix( PDE, Nx, Ny, dx, dy );
%
% Inputs:
%  PDE  PDE object containing information about the system
%  Nx - number of discretization columns on the x axis
%  Ny - number of discretization rows on the y axis
%  dx - step size in the Ox direction
%  dy - step size in the Oy direction
%
% Outputs:
%  R  matrix of size Ny x Nx by Ny x Nx, where Nx is the number of columns
%     and Ny the number of rows of the mesh


% R = sparse(Ny*Nx, Ny*Nx);

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
                K_12 =   0.5 * PDE.K.eval_coef(1,2,x_pos-dx,y_pos) + ...
                       + 0.5 * PDE.K.eval_coef(1,2,x_pos,y_pos);
                K_21 =   0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos-dy) + ...
                       + 0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos);
                delta = K_12*dx*dy + K_21*dx*dy;
                R(k,k-Nx-1) = delta;
            end
            K_12_left =   0.5 * PDE.K.eval_coef(1,2,x_pos-dx,y_pos) + ...
                        + 0.5 * PDE.K.eval_coef(1,2,x_pos,y_pos);
            K_12_right =   0.5 * PDE.K.eval_coef(1,2,x_pos+dx,y_pos) + ...
                         + 0.5 * PDE.K.eval_coef(1,2,x_pos,y_pos);
            K_22 =   0.5 * PDE.K.eval_coef(2,2,x_pos,y_pos-dy) + ...
                   + 0.5 * PDE.K.eval_coef(2,2,x_pos,y_pos);
            gamma = K_12_left*dx*dy - K_12_right*dx*dy + 4*K_22*dx^2;
            R(k,k-Nx) = gamma;
            if nx < Nx
                K_12 =   0.5 * PDE.K.eval_coef(1,2,x_pos+dx,y_pos) + ...
                       + 0.5 * PDE.K.eval_coef(1,2,x_pos,y_pos);
                K_21 =   0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos-dy) + ...
                       + 0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos);
                delta = K_12*dx*dy + K_21*dx*dy;
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
            K_21_below =   0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos-dy) + ...
                         + 0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos);
            K_21_above =   0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos+dy) + ...
                         + 0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos);
            K_11 =   0.5 * PDE.K.eval_coef(1,1,x_pos-dx,y_pos) + ...
                   + 0.5 * PDE.K.eval_coef(1,1,x_pos,y_pos);
            beta =   K_21_below*dx*dy - K_21_above*dx*dy + 4*K_11*dy^2;
            R(k,k-1) = beta;
        end
        K_11_left =   0.5 * PDE.K.eval_coef(1,1,x_pos-dx,y_pos) + ...
                    + 0.5 * PDE.K.eval_coef(1,1,x_pos,y_pos);
        K_11_right =   0.5 * PDE.K.eval_coef(1,1,x_pos+dx,y_pos) + ...
                     + 0.5 * PDE.K.eval_coef(1,1,x_pos,y_pos);
        K_22_below =   0.5 * PDE.K.eval_coef(2,2,x_pos,y_pos-dy) + ...
                     + 0.5 * PDE.K.eval_coef(2,2,x_pos,y_pos);
        K_22_above =   0.5 * PDE.K.eval_coef(2,2,x_pos,y_pos+dy) + ...
                     + 0.5 * PDE.K.eval_coef(2,2,x_pos,y_pos);
        alpha = - 4*K_11_left*dy^2 - 4*K_11_right*dy^2 + ...
                - 4*K_22_below*dx^2 - 4*K_22_above*dx^2;
        R(k,k) = alpha;
        if nx < Nx
            K_21_above =   0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos+dy) + ...
                         + 0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos);
            K_21_below =   0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos-dy) + ...
                         + 0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos);
            K_11 =   0.5 * PDE.K.eval_coef(1,1,x_pos+dx,y_pos) + ...
                   + 0.5 * PDE.K.eval_coef(1,1,x_pos,y_pos);
            beta = K_21_above*dx*dy - K_21_below*dx*dy + 4*K_11*dy^2;
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
                K_12 =   0.5 * PDE.K.eval_coef(1,2,x_pos-dx,y_pos) + ...
                       + 0.5 * PDE.K.eval_coef(1,2,x_pos,y_pos);
                K_21 =   0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos+dy) + ...
                       + 0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos);
                delta = K_12*dx*dy + K_21*dx*dy;
                R(k,k+Nx-1) = -delta;
            end
            K_12_right =   0.5 * PDE.K.eval_coef(1,2,x_pos+dx,y_pos) + ...
                         + 0.5 * PDE.K.eval_coef(1,2,x_pos,y_pos);
            K_12_left =   0.5 * PDE.K.eval_coef(1,2,x_pos-dx,y_pos) + ...
                        + 0.5 * PDE.K.eval_coef(1,2,x_pos,y_pos);
            K_22 =   0.5 * PDE.K.eval_coef(2,2,x_pos,y_pos+dy) + ...
                   + 0.5 * PDE.K.eval_coef(2,2,x_pos,y_pos);
            gamma = K_12_right*dx*dy - K_12_left*dx*dy + 4*K_22*dx^2;
            R(k,k+Nx) = gamma;
            if nx < Nx
                K_12 =   0.5 * PDE.K.eval_coef(1,2,x_pos+dx,y_pos) + ...
                       + 0.5 * PDE.K.eval_coef(1,2,x_pos,y_pos);
                K_21 =   0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos+dy) + ...
                       + 0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos);
                delta = K_12*dx*dy + K_21*dx*dy;
                R(k,k+Nx+1) = delta;
            end
        end
    end

end

end