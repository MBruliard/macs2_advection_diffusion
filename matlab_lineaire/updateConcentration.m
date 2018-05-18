function [ C_updated ] = updateConcentration( PDE, C, V, Nx, Ny, dt )
% updateConcentration: return the concentration matrix C given the PDE
% system and the curent state concentration.
% [ C ] = updateConcentration( PDE, C, V, Nx, Ny, dt );
%
% Inputs:
%  PDE: PDE object containing information about the system
%  C:   concentration matrix of size Ny by Nx, where Nx is the number of
%       columns and Ny the number of rows of the mesh
%  V.x: velocity matrix (Ox direction) of size Nx x Ny by 4, where the
%       first argument is the cell number and the second one is the border
%       indice (anticlockwise, lower border is number 1)
%  V.y: velocity matrix (Oy direction) of size Nx x Ny by 4, where the
%       first argument is the cell number and the second one is the border
%       indice (anticlockwise, lower border is number 1)
%  Nx:  number of discretization columns on the x axis
%  Ny:  number of discretization rows on the y axis
%  dt:  step of time
%
% Outputs:
%  C_updated: concentration matrix of size Ny by Nx, where Nx is the number
%             of columns and Ny the number of rows of the mesh


% Cell borders
LOWER_BORDER = 1;
RIGHT_BORDER = 2;
UPPER_BORDER = 3;
LEFT_BORDER = 4;

dx = (PDE.x_max - PDE.x_min) / (Nx-1);
dy = (PDE.y_max - PDE.y_min) / (Ny-1);

C_updated = zeros(Ny, Nx);

% We update the concentration in the inner part of the mesh
for i = 2 : Nx-1
    for j = 2 : Ny-1
        k = i + (j-1)*Nx;

        if V.y(k, LOWER_BORDER) <= 0
            lower_flux = V.y(k, LOWER_BORDER) * C(j,i);
        else
            lower_flux = V.y(k, LOWER_BORDER) * C(j-1,i);
        end

        if V.x(k, RIGHT_BORDER) >= 0
            right_flux = V.x(k, RIGHT_BORDER) *C(j,i);
        else
            right_flux = V.x(k, RIGHT_BORDER) * C(j,i+1);
        end

        if V.y(k, UPPER_BORDER) >= 0
            upper_flux = V.y(k, UPPER_BORDER) * C(j,i);
        else
            upper_flux = V.y(k, UPPER_BORDER) * C(j+1,i);
        end

        if V.x(k, LEFT_BORDER) <= 0
            left_flux = V.x(k, LEFT_BORDER) * C(j,i);
        else
            left_flux = V.x(k, LEFT_BORDER) * C(j,i-1);
        end

        C_updated(j,i) = C(j,i) + ...
            - dt/ dx * ( right_flux - left_flux) + ...
            - dt/ dy * ( upper_flux - lower_flux);

    end
end


% % % FIRST ORDER BORDERS
% % We update the concentration on the lower and upper borders
% for j = 2 : Nx-1
%     C_updated(1,j)  = C_updated(2,j);
%     C_updated(Ny,j) = C_updated(Ny-1,j);
% end
% 
% % We update the concentration on the left and right borders
% for i = 2 : Ny-1
%     C_updated(i,1) =  C_updated(i,2);
%     C_updated(i,Nx) = C_updated(i,Nx-1);
% end
% 
% % We update the concentration on the lower left corner
% C_updated(1,1) = C_updated(2,2);
% 
% % We update the concentration on the lower right corner
% C_updated(1,Nx) = C_updated(2,Nx-1);
% 
% % We update the concentration on the upper left corner
% C_updated(Ny,1) = C_updated(Ny-1,2);
% 
% % We update the concentration on the upper right corner
% C_updated(Ny,Nx) = C_updated(Ny-1,Nx-1);


% % SECOND ORDER BORDERS
% We update the concentration on the lower and upper borders
for j = 2 : Nx-1
    C_updated(1,j)  = (4 * C_updated(2,j) - C_updated(3,j)) / 3;
    C_updated(Ny,j) = (4 * C_updated(Ny-1,j) - C_updated(Ny-2,j)) / 3;
end

% We update the concentration on the left and right borders
for i = 2 : Ny-1
    C_updated(i,1) = (4 * C_updated(i,2) - C_updated(i,3)) / 3;
    C_updated(i,Nx) = (4 * C_updated(i,Nx-1) - C_updated(i,Nx-2)) / 3;
end

% We update the concentration on the lower left corner
C_updated(1,1) = (4 * C_updated(2,2) - C_updated(3,3)) / 3;

% We update the concentration on the lower right corner
C_updated(1,Nx) = (4 * C_updated(2,Nx-1) - C_updated(3,Nx-2)) / 3;

% We update the concentration on the upper left corner
C_updated(Ny,1) = (4 * C_updated(Ny-1,2) - C_updated(Ny-2,3)) / 3;

% We update the concentration on the upper right corner
C_updated(Ny,Nx) = (4 * C_updated(Ny-1,Nx-1) - C_updated(Ny-2,Nx-2)) / 3;


end