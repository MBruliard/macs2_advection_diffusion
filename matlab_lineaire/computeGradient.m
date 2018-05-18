function [ Grad ] = computeGradient( PDE, H, Nx, Ny )
% computeGradient: return the gradient of the physical quantity H
% [ Grad ] = computeGradient( PDE, H, Nx, Ny );
%
% Inputs:
%  PDE  PDE object containing information about the system
%  H -- field of a physical quantity (matrix of size Ny by Nx)
%  Nx - number of discretization columns on the x axis
%  Ny - number of discretization rows on the y axis
%
% Outputs:
%  Grad -- gradient of the physical quantity with components Grad.x and
%          Grad.y
%  Grad.x  matrix of size Nx x Ny by 4 representing the quantity variation
%          (Ox direction), where the first argument is the cell number
%          and the second one is the border index (anticlockwise, lower
%          border is number 1)
%  Grad.y  matrix of size Nx x Ny by 4 representing the quantity variation
%          (Oy direction), where the first argument is the cell number
%          and the second one is the border index (anticlockwise, lower
%          border is number 1)


% Cell borders
LOWER_BORDER = 1;
RIGHT_BORDER = 2;
UPPER_BORDER = 3;
LEFT_BORDER  = 4;

dx = (PDE.x_max - PDE.x_min) / (Nx-1);
dy = (PDE.y_max - PDE.y_min) / (Ny-1);

Grad.x = zeros(Nx*Ny, 4); % Quantity variation (Ox direction)
Grad.y = zeros(Nx*Ny, 4); % Quantity variation (Oy direction)

% We set the values of GradH.x and GradH.y on each border
for i = 2 : Nx-1
    for j = 2 : Ny-1
        k = i + (j-1)*Nx;

        Grad.x(k, LOWER_BORDER) = (H(j,i+1) + H(j-1,i+1) - H(j,i-1) - H(j-1,i-1)) / (4*dx);
        Grad.x(k, RIGHT_BORDER) = (H(j,i+1) - H(j,i)) / dx;
        Grad.x(k, UPPER_BORDER) = (H(j,i+1) + H(j+1,i+1) - H(j,i-1) - H(j+1,i-1)) / (4*dx);
        Grad.x(k, LEFT_BORDER)  = (H(j,i) - H(j,i-1)) / dx;

        Grad.y(k, LOWER_BORDER) = (H(j,i) - H(j-1,i)) / dy;
        Grad.y(k, RIGHT_BORDER) = (H(j+1,i+1) + H(j+1,i) - H(j-1,i+1) - H(j-1,i)) / (4*dy);
        Grad.y(k, UPPER_BORDER) = (H(j+1,i) - H(j,i)) / dy;
        Grad.y(k, LEFT_BORDER)  = (H(j+1,i-1) + H(j+1,i) - H(j-1,i-1) - H(j-1,i)) / (4*dy);
    end
end

% Lower domain border
for i = 2 : Nx-1
    k = i;

    Grad.x(k, RIGHT_BORDER) = (H(1,i+1) - H(1,i)) / dx;
    Grad.x(k, UPPER_BORDER) = (H(1,i+1) + H(2,i+1) - H(1,i-1) - H(2,i-1)) / (4*dx);
    Grad.x(k, LEFT_BORDER)  = (H(1,i) - H(1,i-1)) / dx;
    Grad.x(k, LOWER_BORDER) = (Grad.x(k, LEFT_BORDER) + Grad.x(k, RIGHT_BORDER)) / 2;

    Grad.y(k, RIGHT_BORDER) = (H(2,i+1) + H(2,i) - H(1,i+1) - H(1,i)) / (2*dy);
    Grad.y(k, UPPER_BORDER) = (H(2,i) - H(1,i)) / dy;
    Grad.y(k, LEFT_BORDER)  = (H(2,i) + H(2,i-1) - H(1,i) - H(1,i-1)) / (2*dy);
    Grad.y(k, LOWER_BORDER) = (Grad.y(k, LEFT_BORDER) + Grad.y(k, RIGHT_BORDER)) / 2;
end

% Right domain border
for j = 2 : Ny-1
    k = Nx + (j-1)*Nx;
    Grad.x(k, LOWER_BORDER) = (H(j,Nx) + H(j-1,Nx) - H(j,Nx-1) - H(j-1,Nx-1)) / (2*dx);
    Grad.x(k, UPPER_BORDER) = (H(j,Nx) + H(j+1,Nx) - H(j,Nx-1) - H(j+1,Nx-1)) / (2*dx);
    Grad.x(k, LEFT_BORDER)  = (H(j,Nx) - H(j,Nx-1)) / dx;
    Grad.x(k, RIGHT_BORDER) = (Grad.x(k, LOWER_BORDER) + Grad.x(k, UPPER_BORDER)) / 2;
    
    Grad.y(k, LOWER_BORDER) = (H(j,Nx) - H(j-1,Nx)) / dy;
    Grad.y(k, UPPER_BORDER) = (H(j+1,Nx) - H(j,Nx)) / dy;
    Grad.y(k, LEFT_BORDER)  = (H(j+1,Nx) + H(j+1,Nx-1) - H(j-1,Nx) - H(j-1,Nx-1)) / (4*dy);
    Grad.y(k, RIGHT_BORDER) = (Grad.y(k, LOWER_BORDER) + Grad.y(k, UPPER_BORDER)) / 2;
end

% Upper domain border
for i = 2 : Nx-1
    k = i + (Ny-1)*Nx;
    Grad.x(k, LOWER_BORDER) = (H(Ny,i+1) + H(Ny-1,i+1) - H(Ny,i-1) - H(Ny-1,i-1)) / (4*dx);
    Grad.x(k, RIGHT_BORDER) = (H(Ny,i+1) - H(Ny,i)) / dx;
    Grad.x(k, LEFT_BORDER)  = (H(Ny,i) - H(Ny,i-1)) / dx;
    Grad.x(k, UPPER_BORDER) = (Grad.x(k, LEFT_BORDER) + Grad.x(k, RIGHT_BORDER)) / 2;
    
    Grad.y(k, LOWER_BORDER) = (H(Ny,i) - H(Ny-1,i)) / dy;
    Grad.y(k, RIGHT_BORDER) = (H(Ny,i+1) + H(Ny,i) - H(Ny-1,i+1) - H(Ny-1,i)) / (2*dy);
    Grad.y(k, LEFT_BORDER)  = (H(Ny,i) + H(Ny,i-1) - H(Ny-1,i) - H(Ny-1,i-1)) / (2*dy);
    Grad.y(k, UPPER_BORDER) = (Grad.y(k, LEFT_BORDER) + Grad.y(k, RIGHT_BORDER)) / 2;
end

% Left domain border
for j = 2 : Ny-1
    k = 1 + (j-1)*Nx;
    Grad.x(k, LOWER_BORDER) = (H(j,2) + H(j-1,2) - H(j,1) - H(j-1,1)) / (2*dx);
    Grad.x(k, RIGHT_BORDER) = (H(j,2) - H(j,1)) / dx;
    Grad.x(k, UPPER_BORDER) = (H(j,2) + H(j+1,2) - H(j,1) - H(j+1,1)) / (2*dx);
    Grad.x(k, LEFT_BORDER)  = (Grad.x(k, LOWER_BORDER) + Grad.x(k, UPPER_BORDER)) / 2;

    Grad.y(k, LOWER_BORDER) = (H(j,1) - H(j-1,1)) / dy;
    Grad.y(k, RIGHT_BORDER) = (H(j+1,2) + H(j+1,1) - H(j-1,2) - H(j-1,1)) / (4*dy);
    Grad.y(k, UPPER_BORDER) = (H(j+1,1) - H(j,1)) / dy;
    Grad.y(k, LEFT_BORDER)  = (Grad.y(k, LOWER_BORDER) + Grad.y(k, UPPER_BORDER)) / 2;
end

% Lower left corner
Grad.x(1, RIGHT_BORDER) = (H(1,2) - H(1,1)) / dx;
Grad.x(1, UPPER_BORDER) = (H(2,2) + H(1,2) - H(2,1) - H(1,1)) / (2*dx);
Grad.x(1, LOWER_BORDER) = Grad.x(1, RIGHT_BORDER);
Grad.x(1, LEFT_BORDER)  = Grad.x(1, UPPER_BORDER);

Grad.y(1, RIGHT_BORDER) = (H(2,2) + H(2,1) - H(1,2) - H(1,1)) / (2*dy);
Grad.y(1, UPPER_BORDER) = (H(2,1) - H(1,1)) / dy;
Grad.y(1, LOWER_BORDER) = Grad.y(1, RIGHT_BORDER);
Grad.y(1, LEFT_BORDER)  = Grad.y(1, RIGHT_BORDER);

% Lower right corner
Grad.x(Nx, LEFT_BORDER)  = (H(1,Nx) - H(1,Nx-1)) / dx;
Grad.x(Nx, UPPER_BORDER) = (H(2,Nx) + H(1,Nx) - H(2,Nx-1) - H(1,Nx-1)) / (2*dx);
Grad.x(Nx, LOWER_BORDER) = Grad.x(Nx, LEFT_BORDER);
Grad.x(Nx, RIGHT_BORDER) = Grad.x(Nx, UPPER_BORDER);

Grad.y(Nx, LEFT_BORDER)  = (H(2,Nx) + H(2,Nx-1) - H(1,Nx) - H(1,Nx-1)) / (2*dy);
Grad.y(Nx, UPPER_BORDER) = (H(2,Nx) - H(1,Nx)) / dy;
Grad.y(Nx, LOWER_BORDER) = Grad.y(Nx, LEFT_BORDER);
Grad.y(Nx, RIGHT_BORDER) = Grad.y(Nx, UPPER_BORDER);

% Upper right corner
Grad.x(Nx*Ny, LEFT_BORDER)  = (H(Ny,Nx) - H(Ny,Nx-1)) / dx;
Grad.x(Nx*Ny, LOWER_BORDER) = (H(Ny,Nx) + H(Ny-1,Nx) - H(Ny,Nx-1) - H(Ny-1,Nx-1)) / (2*dx);
Grad.x(Nx*Ny, UPPER_BORDER) = Grad.x(Nx*Ny, LEFT_BORDER);
Grad.x(Nx*Ny, RIGHT_BORDER) = Grad.x(Nx*Ny, LOWER_BORDER);

Grad.y(Nx*Ny, LOWER_BORDER) = (H(Ny,Nx) - H(Ny-1,Nx)) / dy;
Grad.y(Nx*Ny, LEFT_BORDER)  = (H(Ny,Nx) + H(Ny,Nx-1) - H(Ny-1,Nx) - H(Ny-1,Nx-1)) / (2*dy);
Grad.y(Nx*Ny, RIGHT_BORDER) = Grad.y(Nx*Ny, LOWER_BORDER);
Grad.y(Nx*Ny, UPPER_BORDER) = Grad.y(Nx*Ny, LEFT_BORDER);

% Upper left corner
Grad.x(1+(Ny-1)*Nx, LOWER_BORDER) = (H(Ny,2) + H(Ny-1,2) - H(Ny,1) - H(Ny-1,1)) / (2*dx);
Grad.x(1+(Ny-1)*Nx, RIGHT_BORDER) = (H(Ny,2) - H(Ny,1)) / dx;
Grad.x(1+(Ny-1)*Nx, UPPER_BORDER) = Grad.x(1+(Ny-1)*Nx, RIGHT_BORDER);
Grad.x(1+(Ny-1)*Nx, LEFT_BORDER)  = Grad.x(1+(Ny-1)*Nx, LOWER_BORDER);

Grad.y(1+(Ny-1)*Nx, LOWER_BORDER) = (H(Ny,1) - H(Ny-1,1)) / dy;
Grad.y(1+(Ny-1)*Nx, RIGHT_BORDER) = (H(Ny,2) + H(Ny,1) - H(Ny-1,2) - H(Ny-1,1)) / (2*dy);
Grad.y(1+(Ny-1)*Nx, LEFT_BORDER)  = Grad.y(1+(Ny-1)*Nx, LOWER_BORDER);
Grad.y(1+(Ny-1)*Nx, UPPER_BORDER) = Grad.y(1+(Ny-1)*Nx, RIGHT_BORDER);

end