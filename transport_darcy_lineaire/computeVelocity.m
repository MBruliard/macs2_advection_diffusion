function [ V ] = computeVelocity( PDE, H )
% computeVelocity: return two velocity fields: one in the X component, and
% the other in the Y.
% [ V ] = computeVelocity( PDE, H, Nx, Ny );
%
% Inputs:
%  PDE  PDE object containing information about the system
%  H -- pressure matrix of size Ny by Nx, where Nx is the number of columns
%       and Ny the number of rows of the mesh
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

GradH = computeGradient( PDE, H );

V.x = zeros(PDE.Ny*PDE.Nx,2);
V.y = zeros(PDE.Ny*PDE.Nx,2);

ny_index_ = 2 : PDE.Ny-1;


for nx = 2 : PDE.Nx-1
    
    k_ = nx + (ny_index_-1) * PDE.Nx;
    
    K_11_below_ = PDE.K.eval_harmonic_mean_(1,1,nx,ny_index_,[0, -1]);
    K_12_below_ = PDE.K.eval_harmonic_mean_(1,2,nx,ny_index_,[0, -1]);
    K_21_below_ = PDE.K.eval_harmonic_mean_(2,1,nx,ny_index_,[0, -1]);
    K_22_below_ = PDE.K.eval_harmonic_mean_(2,2,nx,ny_index_,[0, -1]);
    
    V.x(k_, LOWER_BORDER) = - K_11_below_ .* GradH.x(k_, LOWER_BORDER) + ...
        - K_12_below_ .* GradH.y(k_, LOWER_BORDER);
    V.y(k_, LOWER_BORDER) = - K_21_below_ .* GradH.x(k_, LOWER_BORDER) + ...
        - K_22_below_ .* GradH.y(k_, LOWER_BORDER);

    K_11_right_ = PDE.K.eval_harmonic_mean_(1,1,nx,ny_index_,[1, 0]);
    K_12_right_ = PDE.K.eval_harmonic_mean_(1,2,nx,ny_index_,[1, 0]);
    K_21_right_ = PDE.K.eval_harmonic_mean_(2,1,nx,ny_index_,[1, 0]);
    K_22_right_ = PDE.K.eval_harmonic_mean_(2,2,nx,ny_index_,[1, 0]);

    V.x(k_, RIGHT_BORDER) = - K_11_right_ .* GradH.x(k_, RIGHT_BORDER) + ...
        - K_12_right_ .* GradH.y(k_, RIGHT_BORDER);
    V.y(k_, RIGHT_BORDER) = - K_21_right_ .* GradH.x(k_, RIGHT_BORDER) + ...
        - K_22_right_ .* GradH.y(k_, RIGHT_BORDER);

end


% Lower border
nx = 1;
for ny = 2 : PDE.Nx-1
    k = ny + (nx-1)*PDE.Nx;
    V.x(k, LOWER_BORDER) = V.x(k+PDE.Nx, LOWER_BORDER);
    V.y(k, LOWER_BORDER) = V.y(k+PDE.Nx, LOWER_BORDER);
    V.x(k, RIGHT_BORDER) = V.x(k+PDE.Nx, RIGHT_BORDER);
    V.y(k, RIGHT_BORDER) = V.y(k+PDE.Nx, RIGHT_BORDER);
end

% Right border
nx = PDE.Nx;
for ny = 2 : PDE.Ny-1
    k = nx + (ny-1)*PDE.Nx;
    V.x(k, LOWER_BORDER) = V.x(k-1, LOWER_BORDER);
    V.y(k, LOWER_BORDER) = V.y(k-1, LOWER_BORDER);
    V.x(k, RIGHT_BORDER) = V.x(k-1, RIGHT_BORDER);
    V.y(k, RIGHT_BORDER) = V.y(k-1, RIGHT_BORDER);
end

% Upper border
nx = PDE.Ny;
for ny = 2 : PDE.Nx-1
    k = ny + (nx-1)*PDE.Nx;
    V.x(k, LOWER_BORDER) = V.x(k-PDE.Nx, LOWER_BORDER);
    V.y(k, LOWER_BORDER) = V.y(k-PDE.Nx, LOWER_BORDER);
    V.x(k, RIGHT_BORDER) = V.x(k-PDE.Nx, RIGHT_BORDER);
    V.y(k, RIGHT_BORDER) = V.y(k-PDE.Nx, RIGHT_BORDER);
end

% Left border
nx = 1;
for ny = 2 : PDE.Ny-1
    k = nx + (ny-1)*PDE.Nx;
    V.x(k, LOWER_BORDER) = V.x(k+1, LOWER_BORDER);
    V.y(k, LOWER_BORDER) = V.y(k+1, LOWER_BORDER);
    V.x(k, RIGHT_BORDER) = V.x(k+1, RIGHT_BORDER);
    V.y(k, RIGHT_BORDER) = V.y(k+1, RIGHT_BORDER);
end

% Lower left cell
nx = 1;
ny = 1;
k = nx + (ny-1)*PDE.Nx;
V.x(k, LOWER_BORDER) = V.x(k+PDE.Nx+1, LOWER_BORDER);
V.y(k, LOWER_BORDER) = V.y(k+PDE.Nx+1, LOWER_BORDER);
V.x(k, RIGHT_BORDER) = V.x(k+PDE.Nx+1, RIGHT_BORDER);
V.y(k, RIGHT_BORDER) = V.y(k+PDE.Nx+1, RIGHT_BORDER);

% Lower right cell
nx = PDE.Nx;
ny = 1;
k = nx + (ny-1)*PDE.Nx;
V.x(k, LOWER_BORDER) = V.x(k+PDE.Nx-1, LOWER_BORDER);
V.y(k, LOWER_BORDER) = V.y(k+PDE.Nx-1, LOWER_BORDER);
V.x(k, RIGHT_BORDER) = V.x(k+PDE.Nx-1, RIGHT_BORDER);
V.y(k, RIGHT_BORDER) = V.y(k+PDE.Nx-1, RIGHT_BORDER);

% Upper right cell
nx = PDE.Nx;
ny = PDE.Ny;
k = nx + (ny-1)*PDE.Nx;
V.x(k, LOWER_BORDER) = V.x(k-PDE.Nx-1, LOWER_BORDER);
V.y(k, LOWER_BORDER) = V.y(k-PDE.Nx-1, LOWER_BORDER);
V.x(k, RIGHT_BORDER) = V.x(k-PDE.Nx-1, RIGHT_BORDER);
V.y(k, RIGHT_BORDER) = V.y(k-PDE.Nx-1, RIGHT_BORDER);

% Upper left cell
nx = 1;
ny = PDE.Ny;
k = nx + (ny-1)*PDE.Nx;
V.x(k, LOWER_BORDER) = V.x(k-PDE.Nx+1, LOWER_BORDER);
V.y(k, LOWER_BORDER) = V.y(k-PDE.Nx+1, LOWER_BORDER);
V.x(k, RIGHT_BORDER) = V.x(k-PDE.Nx+1, RIGHT_BORDER);
V.y(k, RIGHT_BORDER) = V.y(k-PDE.Nx+1, RIGHT_BORDER);


end