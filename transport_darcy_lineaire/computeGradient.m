function [ Grad ] = computeGradient( PDE, H )
% computeGradient: return the gradient of the physical quantity H
% [ Grad ] = computeGradient( PDE, H, Nx, Ny );
%
% Inputs:
%  PDE  PDE object containing information about the system
%  H -- field of a physical quantity (matrix of size Ny by Nx)
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

Nx = PDE.Nx;
Ny = PDE.Ny;

dx = (PDE.x_max - PDE.x_min) / (Nx-1);
dy = (PDE.y_max - PDE.y_min) / (Ny-1);

Grad.x = zeros(Nx*Ny, 2); % Quantity variation (Ox direction)
Grad.y = zeros(Nx*Ny, 2); % Quantity variation (Oy direction)

% We set the values of Grad.x and Grad.y on each border
ny_ = 2 : Ny-1;
for nx = 2 : Nx-1
    k_ = nx + (ny_-1)*Nx;
    Grad.x(k_, LOWER_BORDER) = (H(ny_,nx+1) + H(ny_-1,nx+1) - H(ny_,nx-1) - H(ny_-1,nx-1)) / (4*dx);
    Grad.x(k_, RIGHT_BORDER) = (H(ny_,nx+1) - H(ny_,nx)) / dx;
    Grad.y(k_, LOWER_BORDER) = (H(ny_,nx) - H(ny_-1,nx)) / dy;
    Grad.y(k_, RIGHT_BORDER) = (H(ny_+1,nx+1) + H(ny_+1,nx) - H(ny_-1,nx+1) - H(ny_-1,nx)) / (4*dy);
end

% Lower domain border
for nx = 2 : Nx-1
    k = nx;
    Grad.x(k, LOWER_BORDER) = (Grad.x(k-1, RIGHT_BORDER) + Grad.x(k, RIGHT_BORDER)) / 2;
    Grad.x(k, RIGHT_BORDER) = (H(1,nx+1) - H(1,nx)) / dx;
    Grad.y(k, LOWER_BORDER) = (Grad.y(k-1, RIGHT_BORDER) + Grad.y(k, RIGHT_BORDER)) / 2;
    Grad.y(k, RIGHT_BORDER) = (H(2,nx+1) + H(2,nx) - H(1,nx+1) - H(1,nx)) / (2*dy);
end

% Right domain border
for ny = 2 : Ny-1
    k = Nx + (ny-1)*Nx;
    Grad.x(k, LOWER_BORDER) = (H(ny,Nx) + H(ny-1,Nx) - H(ny,Nx-1) - H(ny-1,Nx-1)) / (2*dx);
    Grad.x(k, RIGHT_BORDER) = (Grad.x(k, LOWER_BORDER) + Grad.x(k+Nx, LOWER_BORDER)) / 2;
    Grad.y(k, LOWER_BORDER) = (H(ny,Nx) - H(ny-1,Nx)) / dy;
    Grad.y(k, RIGHT_BORDER) = (Grad.y(k, LOWER_BORDER) + Grad.y(k+Nx, LOWER_BORDER)) / 2;
end

% Upper domain border
for nx = 2 : Nx-1
    k = nx + (Ny-1)*Nx;
    Grad.x(k, LOWER_BORDER) = (H(Ny,nx+1) + H(Ny-1,nx+1) - H(Ny,nx-1) - H(Ny-1,nx-1)) / (4*dx);
    Grad.x(k, RIGHT_BORDER) = (H(Ny,nx+1) - H(Ny,nx)) / dx;
    Grad.y(k, LOWER_BORDER) = (H(Ny,nx) - H(Ny-1,nx)) / dy;
    Grad.y(k, RIGHT_BORDER) = (H(Ny,nx+1) + H(Ny,nx) - H(Ny-1,nx+1) - H(Ny-1,nx)) / (2*dy);
end

% Left domain border
for ny = 2 : Ny-1
    k = 1 + (ny-1)*Nx;
    Grad.x(k, LOWER_BORDER) = (H(ny,2) + H(ny-1,2) - H(ny,1) - H(ny-1,1)) / (2*dx);
    Grad.x(k, RIGHT_BORDER) = (H(ny,2) - H(ny,1)) / dx;
    Grad.y(k, LOWER_BORDER) = (H(ny,1) - H(ny-1,1)) / dy;
    Grad.y(k, RIGHT_BORDER) = (H(ny+1,2) + H(ny+1,1) - H(ny-1,2) - H(ny-1,1)) / (4*dy);
end

% Lower left corner
Grad.x(1, LOWER_BORDER) = Grad.x(1, RIGHT_BORDER);
Grad.x(1, RIGHT_BORDER) = (H(1,2) - H(1,1)) / dx;
Grad.y(1, LOWER_BORDER) = Grad.y(1, RIGHT_BORDER);
Grad.y(1, RIGHT_BORDER) = (H(2,2) + H(2,1) - H(1,2) - H(1,1)) / (2*dy);

% Lower right corner
Grad.x(Nx, LOWER_BORDER) = Grad.x(Nx-1, RIGHT_BORDER);
Grad.x(Nx, RIGHT_BORDER) = Grad.x(2*Nx, LOWER_BORDER);
Grad.y(Nx, LOWER_BORDER) = Grad.y(Nx-1, RIGHT_BORDER);
Grad.y(Nx, RIGHT_BORDER) = Grad.y(2*Nx, LOWER_BORDER);

% Upper right corner
Grad.x(Nx*Ny, LOWER_BORDER) = (H(Ny,Nx) + H(Ny-1,Nx) - H(Ny,Nx-1) - H(Ny-1,Nx-1)) / (2*dx);
Grad.x(Nx*Ny, RIGHT_BORDER) = Grad.x(Nx*Ny, LOWER_BORDER);
Grad.y(Nx*Ny, LOWER_BORDER) = (H(Ny,Nx) - H(Ny-1,Nx)) / dy;
Grad.y(Nx*Ny, RIGHT_BORDER) = Grad.y(Nx*Ny, LOWER_BORDER);

% Upper left corner
Grad.x(1+(Ny-1)*Nx, LOWER_BORDER) = (H(Ny,2) + H(Ny-1,2) - H(Ny,1) - H(Ny-1,1)) / (2*dx);
Grad.x(1+(Ny-1)*Nx, RIGHT_BORDER) = (H(Ny,2) - H(Ny,1)) / dx;
Grad.y(1+(Ny-1)*Nx, LOWER_BORDER) = (H(Ny,1) - H(Ny-1,1)) / dy;
Grad.y(1+(Ny-1)*Nx, RIGHT_BORDER) = (H(Ny,2) + H(Ny,1) - H(Ny-1,2) - H(Ny-1,1)) / (2*dy);

end