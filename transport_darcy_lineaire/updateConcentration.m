function [ C_updated ] = updateConcentration( PDE, C, V, dt )
% updateConcentration: return the concentration matrix C given the PDE
% system and the curent state concentration.
% [ C ] = updateConcentration( PDE, C, V, Nx, Ny, dt );
%
% Inputs:
%  PDE  PDE object containing information about the system
%  C -- concentration matrix of size Ny by Nx, where Nx is the number
%       of columns and Ny the number of rows of the mesh
%  V.x  velocity matrix (Ox direction) of size Nx x Ny by 4, where the
%       first argument is the cell number and the second one is the
%       border index (anticlockwise, lower border is number 1)
%  V.y  velocity matrix (Oy direction) of size Nx x Ny by 4, where the
%       first argument is the cell number and the second one is the
%       border index (anticlockwise, lower border is number 1)
%  dt - step of time
%
% Outputs:
%  C_updated  concentration matrix of size Ny by Nx, where Nx is the number
%             of columns and Ny the number of rows of the mesh


% Cell borders
LOWER_BORDER = 1;
RIGHT_BORDER = 2;

Nx = PDE.Nx;
Ny = PDE.Ny;

dx = (PDE.x_max - PDE.x_min) / (Nx-1);
dy = (PDE.y_max - PDE.y_min) / (Ny-1);

C_updated = zeros(Ny, Nx);


% We update the concentration in the inner part of the mesh
ny_ = 2 : Ny-1;
for nx = 2 : Nx-1
    k_ = nx + (ny_-1)*Nx;
    lower_flux_ = ...
        V.y(k_, LOWER_BORDER) .* C(ny_,nx) .* (V.y(k_, LOWER_BORDER) <= 0) ...
        + V.y(k_, LOWER_BORDER) .* C(ny_-1,nx) .* (V.y(k_, LOWER_BORDER) > 0);
    right_flux_ = ...
        V.x(k_, RIGHT_BORDER) .* C(ny_,nx) .* (V.x(k_, RIGHT_BORDER) >= 0) ...
        + V.x(k_, RIGHT_BORDER) .* C(ny_,nx+1) .* (V.x(k_, RIGHT_BORDER) < 0);
    upper_flux_ = ...
        V.y(k_+Nx, LOWER_BORDER) .* C(ny_,nx) .* (V.y(k_+Nx, LOWER_BORDER) >= 0) ...
        + V.y(k_+Nx, LOWER_BORDER) .* C(ny_+1,nx) .* (V.y(k_+Nx, LOWER_BORDER) < 0);
    left_flux_ = ...
        V.x(k_-1, RIGHT_BORDER) .* C(ny_,nx) .* (V.x(k_-1, RIGHT_BORDER) <= 0) ...
        + V.x(k_-1, RIGHT_BORDER) .* C(ny_,nx-1) .* (V.x(k_-1, RIGHT_BORDER) > 0);
    C_updated(ny_,nx) = C(ny_,nx) + ...
            - dt / dx * ( right_flux_ - left_flux_) + ...
            - dt / dy * ( upper_flux_ - lower_flux_);
end

if true == PDE.compute_diffusion
    GradC = computeGradient( PDE, C );

    for nx = 2 : Nx-1
        k_ = nx + (ny_-1)*Nx;
        % Computing the diffusion matrix at interfaces
        D_21_below_ = PDE.D.eval_harmonic_mean_(2,1,nx,ny_, [0, -1]);
        D_22_below_ = PDE.D.eval_harmonic_mean_(2,2,nx,ny_, [0, -1]);
        D_11_right_ = PDE.D.eval_harmonic_mean_(1,1,nx,ny_, [1, 0]);
        D_12_right_ = PDE.D.eval_harmonic_mean_(1,2,nx,ny_, [1, 0]);
        D_21_above_ = PDE.D.eval_harmonic_mean_(2,1,nx,ny_, [0, 1]);
        D_22_above_ = PDE.D.eval_harmonic_mean_(2,2,nx,ny_, [0, 1]);
        D_11_left_  = PDE.D.eval_harmonic_mean_(1,1,nx,ny_, [-1, 0]);
        D_12_left_  = PDE.D.eval_harmonic_mean_(1,2,nx,ny_, [-1, 0]);

        % Computing the concentration diffusion through each cell
        lower_diffusion_ = (D_21_below_ .* GradC.x(k_, LOWER_BORDER) + ...
            + D_22_below_ .* GradC.y(k_, LOWER_BORDER)) / dy;
        right_diffusion_ = (D_11_right_ .* GradC.x(k_, RIGHT_BORDER) + ...
            + D_12_right_ .* GradC.y(k_, RIGHT_BORDER)) / dx;
        upper_diffusion_ = (D_21_above_ .* GradC.x(k_+Nx, LOWER_BORDER) + ...
            + D_22_above_ .* GradC.y(k_+Nx, LOWER_BORDER)) / dy;
        left_diffusion_ = (D_11_left_ .* GradC.x(k_-1, RIGHT_BORDER) + ...
            + D_12_left_ .* GradC.y(k_-1, RIGHT_BORDER)) / dx;

        diffusion_ = - lower_diffusion_ + right_diffusion_ + ...
            + upper_diffusion_ - left_diffusion_;
        
        C_updated(ny_,nx) = C_updated(ny_,nx) + dt * diffusion_;
        
    end
end


% Second order borders

% We update the concentration on the lower and upper borders
for ny = 2 : Nx-1
    C_updated(1,ny)  = (4 * C_updated(2,ny) - C_updated(3,ny)) / 3;
    C_updated(Ny,ny) = (4 * C_updated(Ny-1,ny) - C_updated(Ny-2,ny)) / 3;
end

% We update the concentration on the left and right borders
for nx = 2 : Ny-1
    C_updated(nx,1) = (4 * C_updated(nx,2) - C_updated(nx,3)) / 3;
    C_updated(nx,Nx) = (4 * C_updated(nx,Nx-1) - C_updated(nx,Nx-2)) / 3;
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