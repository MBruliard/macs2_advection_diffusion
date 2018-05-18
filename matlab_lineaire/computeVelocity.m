function [ V ] = computeVelocity( PDE, H, Nx, Ny )
% computeVelocity: return two velocity fields: one in the X component, and
% the other in the Y.
% [ V ] = computeVelocity( PDE, H, Nx, Ny );
%
% Inputs:
%  PDE  PDE object containing information about the system
%  H -- pressure matrix of size Ny by Nx, where Nx is the number of columns
%       and Ny the number of rows of the mesh
%  Nx - number of discretization columns on the x axis
%  Ny - number of discretization rows on the y axis
%
% Outputs:
%  V -- velocity object with two subfields V.x and V.y
%  V.x  velocity matrix (Ox direction) of size Nx x Ny by 4, where the
%       first argument is the cell number and the second one is the border
%       indice (anticlockwise, lower border is number 1)
%  V.y  velocity matrix (Oy direction) of size Nx x Ny by 4, where the
%       first argument is the cell number and the second one is the border
%       indice (anticlockwise, lower border is number 1)


% Cell borders
LOWER_BORDER = 1;
RIGHT_BORDER = 2;
UPPER_BORDER = 3;
LEFT_BORDER  = 4;

dx = (PDE.x_max - PDE.x_min) / (Nx-1);
dy = (PDE.y_max - PDE.y_min) / (Ny-1);

GradH = computeGradient( PDE, H, Nx, Ny );

% Darcy's velocity
for i = 1 : Nx
    for j = 1 : Ny
        k = i + (j-1)*Nx;
        x_pos = PDE.x_min + i*dx; % Current x position
        y_pos = PDE.y_min + j*dy; % Current y position

        K_11_below =   0.5 * PDE.K.eval_coef(1,1,x_pos,y_pos-dy) + ...
                     + 0.5 * PDE.K.eval_coef(1,1,x_pos,y_pos);
        K_12_below =   0.5 * PDE.K.eval_coef(1,2,x_pos,y_pos-dy) + ...
                     + 0.5 * PDE.K.eval_coef(1,2,x_pos,y_pos);
        K_21_below =   0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos-dy) + ...
                     + 0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos);
        K_22_below =   0.5 * PDE.K.eval_coef(2,2,x_pos,y_pos-dy) + ...
                     + 0.5 * PDE.K.eval_coef(2,2,x_pos,y_pos);

        V.x(k, LOWER_BORDER) = - K_11_below * GradH.x(k, LOWER_BORDER) + ...
                               - K_12_below * GradH.y(k, LOWER_BORDER);
        V.y(k, LOWER_BORDER) = - K_21_below * GradH.x(k, LOWER_BORDER) + ...
                               - K_22_below * GradH.y(k, LOWER_BORDER);
                           
        K_11_right =   0.5 * PDE.K.eval_coef(1,1,x_pos+dx,y_pos) + ...
                     + 0.5 * PDE.K.eval_coef(1,1,x_pos,y_pos);
        K_12_right =   0.5 * PDE.K.eval_coef(1,2,x_pos+dx,y_pos) + ...
                     + 0.5 * PDE.K.eval_coef(1,2,x_pos,y_pos);
        K_21_right =   0.5 * PDE.K.eval_coef(2,1,x_pos+dx,y_pos) + ...
                     + 0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos);
        K_22_right =   0.5 * PDE.K.eval_coef(2,2,x_pos+dx,y_pos) + ...
                     + 0.5 * PDE.K.eval_coef(2,2,x_pos,y_pos);

        V.x(k, RIGHT_BORDER) = - K_11_right * GradH.x(k, RIGHT_BORDER) + ...
                               - K_12_right * GradH.y(k, RIGHT_BORDER);
        V.y(k, RIGHT_BORDER) = - K_21_right * GradH.x(k, RIGHT_BORDER) + ...
                               - K_22_right * GradH.y(k, RIGHT_BORDER);
                           
        K_11_above =   0.5 * PDE.K.eval_coef(1,1,x_pos,y_pos+dy) + ...
                     + 0.5 * PDE.K.eval_coef(1,1,x_pos,y_pos);
        K_12_above =   0.5 * PDE.K.eval_coef(1,2,x_pos,y_pos+dy) + ...
                     + 0.5 * PDE.K.eval_coef(1,2,x_pos,y_pos);
        K_21_above =   0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos+dy) + ...
                     + 0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos);
        K_22_above =   0.5 * PDE.K.eval_coef(2,2,x_pos,y_pos+dy) + ...
                     + 0.5 * PDE.K.eval_coef(2,2,x_pos,y_pos);

        V.x(k, UPPER_BORDER) = - K_11_above * GradH.x(k, UPPER_BORDER) + ...
                               - K_12_above * GradH.y(k, UPPER_BORDER);
        V.y(k, UPPER_BORDER) = - K_21_above * GradH.x(k, UPPER_BORDER) + ...
                               - K_22_above * GradH.y(k, UPPER_BORDER);
                           
        K_11_left =   0.5 * PDE.K.eval_coef(1,1,x_pos-dx,y_pos) + ...
                    + 0.5 * PDE.K.eval_coef(1,1,x_pos,y_pos);
        K_12_left =   0.5 * PDE.K.eval_coef(1,2,x_pos-dx,y_pos) + ...
                    + 0.5 * PDE.K.eval_coef(1,2,x_pos,y_pos);
        K_21_left =   0.5 * PDE.K.eval_coef(2,1,x_pos-dx,y_pos) + ...
                    + 0.5 * PDE.K.eval_coef(2,1,x_pos,y_pos);
        K_22_left =   0.5 * PDE.K.eval_coef(2,2,x_pos-dx,y_pos) + ...
                    + 0.5 * PDE.K.eval_coef(2,2,x_pos,y_pos);

        V.x(k, LEFT_BORDER) = - K_11_left * GradH.x(k, LEFT_BORDER) + ...
                              - K_12_left * GradH.y(k, LEFT_BORDER);
        V.y(k, LEFT_BORDER) = - K_21_left * GradH.x(k, LEFT_BORDER) + ...
                              - K_22_left * GradH.y(k, LEFT_BORDER);

    end
end

end