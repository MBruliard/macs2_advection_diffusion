function [ V ] = computeVelocity( PDE, H, Nx, Ny )
% computeVelocity: return two velocity fields: one in the X component, and
% the other in the Y.
% [ V ] = computeVelocity( PDE, H, Nx, Ny );
%
% Inputs:
%  PDE: PDE object containing information about the system
%  H:   pressure field (matrix of size Ny by Nx)
%  Nx:  number of discretization columns on the x axis
%  Ny:  number of discretization rows on the y axis
%
% Outputs:
%  V:   velocity object with two subfields V.x and V.y
%  V.x: velocity matrix (Ox direction) of size Nx x Ny by 4, where the
%       first argument is the cell number and the second one is the border
%       indice (anticlockwise, lower border is number 1)
%  V.y: velocity matrix (Oy direction) of size Nx x Ny by 4, where the
%       first argument is the cell number and the second one is the border
%       indice (anticlockwise, lower border is number 1)


% Cell borders
LOWER_BORDER = 1;
RIGHT_BORDER = 2;
UPPER_BORDER = 3;
LEFT_BORDER  = 4;

dx = (PDE.x_max - PDE.x_min) / (Nx-1);
dy = (PDE.y_max - PDE.y_min) / (Ny-1);

vx = zeros(Nx*Ny, 4); % Flux velocity (Ox direction)
vy = zeros(Nx*Ny, 4); % Flux velocity (Oy direction)

% We set the values of vx and vy on each border
for i = 2 : Nx-1
    for j = 2 : Ny-1
        k = i + (j-1)*Nx;

        vx(k, LOWER_BORDER) = (H(j,i+1) + H(j-1,i+1) - H(j,i-1) - H(j-1,i-1)) / (4*dx);
        vx(k, RIGHT_BORDER) = (H(j,i+1) - H(j,i)) / dx;
        vx(k, UPPER_BORDER) = (H(j,i+1) + H(j+1,i+1) - H(j,i-1) - H(j+1,i-1)) / (4*dx);
        vx(k, LEFT_BORDER) = (H(j,i) - H(j,i-1)) / dx;

        vy(k, LOWER_BORDER) = (H(j,i) - H(j-1,i)) / dy;
        vy(k, RIGHT_BORDER) = (H(j+1,i+1) + H(j+1,i) - H(j-1,i+1) - H(j-1,i)) / (4*dy);
        vy(k, UPPER_BORDER) = (H(j+1,i) - H(j,i)) / dy;
        vy(k, LEFT_BORDER) = (H(j+1,i-1) + H(j+1,i) - H(j-1,i-1) - H(j-1,i)) / (4*dy);
    end
end

% Lower domain border
for i = 2 : Nx-1
    k = i;

    vx(k, RIGHT_BORDER) = (H(1,i+1) - H(1,i)) / dx;
    vx(k, UPPER_BORDER) = (H(1,i+1) + H(2,i+1) - H(1,i-1) - H(2,i-1)) / (4*dx);
    vx(k, LEFT_BORDER) = (H(1,i) - H(1,i-1)) / dx;
    vx(k, LOWER_BORDER) = (vx(k, LEFT_BORDER) + vx(k, RIGHT_BORDER)) / 2;

    vy(k, RIGHT_BORDER) = (H(2,i+1) + H(2,i) - H(1,i+1) - H(1,i)) / (2*dy);
    vy(k, UPPER_BORDER) = (H(2,i) - H(1,i)) / dy;
    vy(k, LEFT_BORDER) = (H(2,i) + H(2,i-1) - H(1,i) - H(1,i-1)) / (2*dy);
    vy(k, LOWER_BORDER) = (vy(k, LEFT_BORDER) + vy(k, RIGHT_BORDER)) / 2;
end

% Right domain border
for j = 2 : Ny-1
    k = Nx + (j-1)*Nx;
    vx(k, LOWER_BORDER) = (H(j,Nx) + H(j-1,Nx) - H(j,Nx-1) - H(j-1,Nx-1)) / (2*dx);
    vx(k, UPPER_BORDER) = (H(j,Nx) + H(j+1,Nx) - H(j,Nx-1) - H(j+1,Nx-1)) / (2*dx);
    vx(k, LEFT_BORDER) = (H(j,Nx) - H(j,Nx-1)) / dx;
    vx(k, RIGHT_BORDER) = (vx(k, LOWER_BORDER) + vx(k, UPPER_BORDER)) / 2;
    
    vy(k, LOWER_BORDER) = (H(j,Nx) - H(j-1,Nx)) / dy;
    vy(k, UPPER_BORDER) = (H(j+1,Nx) - H(j,Nx)) / dy;
    vy(k, LEFT_BORDER) = (H(j+1,Nx) + H(j+1,Nx-1) - H(j-1,Nx) - H(j-1,Nx-1)) / (4*dy);
    vy(k, RIGHT_BORDER) = (vy(k, LOWER_BORDER) + vy(k, UPPER_BORDER)) / 2;
end

% Upper domain border
for i = 2 : Nx-1
    k = i + (Ny-1)*Nx;
    vx(k, LOWER_BORDER) = (H(Ny,i+1) + H(Ny-1,i+1) - H(Ny,i-1) - H(Ny-1,i-1)) / (4*dx);
    vx(k, RIGHT_BORDER) = (H(Ny,i+1) - H(Ny,i)) / dx;
    vx(k, LEFT_BORDER) = (H(Ny,i) - H(Ny,i-1)) / dx;
    vx(k, UPPER_BORDER) = (vx(k, LEFT_BORDER) + vx(k, RIGHT_BORDER)) / 2;
    
    vy(k, LOWER_BORDER) = (H(Ny,i) - H(Ny-1,i)) / dy;
    vy(k, RIGHT_BORDER) = (H(Ny,i+1) + H(Ny,i) - H(Ny-1,i+1) - H(Ny-1,i)) / (2*dy);
    vy(k, LEFT_BORDER) = (H(Ny,i) + H(Ny,i-1) - H(Ny-1,i) - H(Ny-1,i-1)) / (2*dy);
    vy(k, UPPER_BORDER) = (vy(k, LEFT_BORDER) + vy(k, RIGHT_BORDER)) / 2;
end

% Left domain border
for j = 2 : Ny-1
    k = 1 + (j-1)*Nx;
    vx(k, LOWER_BORDER) = (H(j,2) + H(j-1,2) - H(j,1) - H(j-1,1)) / (2*dx);
    vx(k, RIGHT_BORDER) = (H(j,2) - H(j,1)) / dx;
    vx(k, UPPER_BORDER) = (H(j,2) + H(j+1,2) - H(j,1) - H(j+1,1)) / (2*dx);
    vx(k, LEFT_BORDER) = (vx(k, LOWER_BORDER) + vx(k, UPPER_BORDER)) / 2;

    vy(k, LOWER_BORDER) = (H(j,1) - H(j-1,1)) / dy;
    vy(k, RIGHT_BORDER) = (H(j+1,2) + H(j+1,1) - H(j-1,2) - H(j-1,1)) / (4*dy);
    vy(k, UPPER_BORDER) = (H(j+1,1) - H(j,1)) / dy;
    vy(k, LEFT_BORDER) = (vy(k, LOWER_BORDER) + vy(k, UPPER_BORDER)) / 2;
end

% Lower left corner
vx(1, RIGHT_BORDER) = (H(1,2) - H(1,1)) / dx;
vx(1, UPPER_BORDER) = (H(2,2) + H(1,2) - H(2,1) - H(1,1)) / (2*dx);
vx(1, LOWER_BORDER) = vx(1, RIGHT_BORDER);
vx(1, LEFT_BORDER) = vx(1, UPPER_BORDER);

vy(1, RIGHT_BORDER) = (H(2,2) + H(2,1) - H(1,2) - H(1,1)) / (2*dy);
vy(1, UPPER_BORDER) = (H(2,1) - H(1,1)) / dy;
vy(1, LOWER_BORDER) = vy(1, RIGHT_BORDER);
vy(1, LEFT_BORDER) = vy(1, RIGHT_BORDER);

% Lower right corner
vx(Nx, LEFT_BORDER) = (H(1,Nx) - H(1,Nx-1)) / dx;
vx(Nx, UPPER_BORDER) = (H(2,Nx) + H(1,Nx) - H(2,Nx-1) - H(1,Nx-1)) / (2*dx);
vx(Nx, LOWER_BORDER) = vx(Nx, LEFT_BORDER);
vx(Nx, RIGHT_BORDER) = vx(Nx, UPPER_BORDER);

vy(Nx, LEFT_BORDER) = (H(2,Nx) + H(2,Nx-1) - H(1,Nx) - H(1,Nx-1)) / (2*dy);
vy(Nx, UPPER_BORDER) = (H(2,Nx) - H(1,Nx)) / dy;
vy(Nx, LOWER_BORDER) = vy(Nx, LEFT_BORDER);
vy(Nx, RIGHT_BORDER) = vy(Nx, UPPER_BORDER);

% Upper right corner
vx(Nx*Ny, LEFT_BORDER) = (H(Ny,Nx) - H(Ny,Nx-1)) / dx;
vx(Nx*Ny, LOWER_BORDER) = (H(Ny,Nx) + H(Ny-1,Nx) - H(Ny,Nx-1) - H(Ny-1,Nx-1)) / (2*dx);
vx(Nx*Ny, UPPER_BORDER) = vx(Nx*Ny, LEFT_BORDER);
vx(Nx*Ny, RIGHT_BORDER) = vx(Nx*Ny, LOWER_BORDER);

vy(Nx*Ny, LOWER_BORDER) = (H(Ny,Nx) - H(Ny-1,Nx)) / dy;
vy(Nx*Ny, LEFT_BORDER) = (H(Ny,Nx) + H(Ny,Nx-1) - H(Ny-1,Nx) - H(Ny-1,Nx-1)) / (2*dy);
vy(Nx*Ny, RIGHT_BORDER) = vy(Nx*Ny, LOWER_BORDER);
vy(Nx*Ny, UPPER_BORDER) = vy(Nx*Ny, LEFT_BORDER);

% Upper left corner
vx(1+(Ny-1)*Nx, LOWER_BORDER) = (H(Ny,2) + H(Ny-1,2) - H(Ny,1) - H(Ny-1,1)) / (2*dx);
vx(1+(Ny-1)*Nx, RIGHT_BORDER) = (H(Ny,2) - H(Ny,1)) / dx;
vx(1+(Ny-1)*Nx, UPPER_BORDER) = vx(1+(Ny-1)*Nx, RIGHT_BORDER);
vx(1+(Ny-1)*Nx, LEFT_BORDER) = vx(1+(Ny-1)*Nx, LOWER_BORDER);

vy(1+(Ny-1)*Nx, LOWER_BORDER) = (H(Ny,1) - H(Ny-1,1)) / dy;
vy(1+(Ny-1)*Nx, RIGHT_BORDER) = (H(Ny,2) + H(Ny,1) - H(Ny-1,2) - H(Ny-1,1)) / (2*dy);
vy(1+(Ny-1)*Nx, LEFT_BORDER) = vy(1+(Ny-1)*Nx, LOWER_BORDER);
vy(1+(Ny-1)*Nx, UPPER_BORDER) = vy(1+(Ny-1)*Nx, RIGHT_BORDER);

% Darcy's velocity
for i = 1 : Nx
    for j = 1 : Ny
        k = i + (j-1)*Nx;
        x_pos = PDE.x_min + i*dx; % Current x position
        y_pos = PDE.y_min + j*dy; % Current y position

        V.x(k, LOWER_BORDER) = ...
            - PDE.K.eval_coef(1,1,x_pos,y_pos) * vx(k, LOWER_BORDER) + ...
            - PDE.K.eval_coef(1,2,x_pos,y_pos) * vy(k, LOWER_BORDER);
        V.y(k, LOWER_BORDER) = ... 
            - PDE.K.eval_coef(2,1,x_pos,y_pos) * vx(k, LOWER_BORDER) + ...
            - PDE.K.eval_coef(2,2,x_pos,y_pos) * vy(k, LOWER_BORDER);

        V.x(k, RIGHT_BORDER) = ...
            - PDE.K.eval_coef(1,1,x_pos,y_pos) * vx(k, RIGHT_BORDER) + ...
            - PDE.K.eval_coef(1,2,x_pos,y_pos) * vy(k, RIGHT_BORDER);
        V.y(k, RIGHT_BORDER) = ...
            - PDE.K.eval_coef(2,1,x_pos,y_pos) * vx(k, RIGHT_BORDER) + ...
            - PDE.K.eval_coef(2,2,x_pos,y_pos) * vy(k, RIGHT_BORDER);

        V.x(k, UPPER_BORDER) = ...
            - PDE.K.eval_coef(1,1,x_pos,y_pos) * vx(k, UPPER_BORDER) + ...
            - PDE.K.eval_coef(1,2,x_pos,y_pos) * vy(k, UPPER_BORDER);
        V.y(k, UPPER_BORDER) = ...
            - PDE.K.eval_coef(2,1,x_pos,y_pos) * vx(k, UPPER_BORDER) + ...
            - PDE.K.eval_coef(2,2,x_pos,y_pos) * vy(k, UPPER_BORDER);

        V.x(k, LEFT_BORDER) = ...
            - PDE.K.eval_coef(1,1,x_pos,y_pos) * vx(k, LEFT_BORDER) + ...
            - PDE.K.eval_coef(1,2,x_pos,y_pos) * vy(k, LEFT_BORDER);
        V.y(k, LEFT_BORDER) = ...
            - PDE.K.eval_coef(2,1,x_pos,y_pos) * vx(k, LEFT_BORDER) + ...
            - PDE.K.eval_coef(2,2,x_pos,y_pos) * vy(k, LEFT_BORDER);

%         K = PDE.K.eval_matrix(x_pos-dx, y_pos-dy);
%         V.x(k, :) = -K(1,1) * vx(k, :) - K(1,2) * vy(k, :);
%         V.y(k, :) = -K(2,1) * vx(k, :) - K(2,2) * vy(k, :);
    end
end

end