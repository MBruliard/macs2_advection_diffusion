function [ ] = solveConcentration( PDE, V, Nx, Ny, CFL )
% solveConcentration: return the concentration matrix C given the PDE
% system.
% [ ] = solveConcentration( PDE, V, Nx, Ny );
%
% Inputs:
%  PDE  PDE object containing information about the system
%  Nx - number of discretization columns on the x axis
%  Ny - number of discretization rows on the y axis
%  V -- velocity object with two subfields V.x and V.y
%  V.x  velocity matrix (Ox direction) of size Nx x Ny by 4, where the
%       first argument is the cell number and the second one is the
%       border index (anticlockwise, lower border is number 1)
%  V.y  velocity matrix (Oy direction) of size Nx x Ny by 4, where the
%       first argument is the cell number and the second one is the
%       border index (anticlockwise, lower border is number 1)


dx = (PDE.x_max - PDE.x_min) / (Nx-1);
dy = (PDE.y_max - PDE.y_min) / (Ny-1);

% Time step (CFL consition)
Vx_max = max(max(abs(V.x)));
Vy_max = max(max(abs(V.y)));
dt_advection = dx*dy / (dy*Vx_max + dx*Vy_max);

D = PDE.D.eval_matrix(0, 0);
if (D(1,1) ~= 0 || D(1,2) ~= 0 || D(2,1) ~= 0 || D(2,2) ~= 0) && true == PDE.compute_diffusion
    dt_diffusion = dx^2*dy^2 / (2*D(1,1)*dy^2 + 2*D(2,2)*dx^2 + (D(1,2)+D(2,1))*dx*dy/2);
    % dt_diffusion = 0.5 * dx^2*dy^2 / ((dx^2+dy^2) * max(max(abs(D))));
else
    dt_diffusion = dt_advection;
end

dt = CFL * 0.5 * min(dt_advection, dt_diffusion);

% Initial conditions
% [X, Y] = meshgrid(0 : dx : PDE.a, 0 : dy : PDE.b);
[X, Y] = PDE.generateMesh(Nx, Ny);
C = PDE.initial_condition(X, Y);
% Z_scale = [min(min(C)) max(max(C))];

hold on
% We loop until final time is reached
time = 0;
while time < PDE.T_final
    dt = min(dt, PDE.T_final-time);
    time = time + dt;
    C = updateConcentration( PDE, C, V, Nx, Ny, dt );
    figure(2);
    surf(X, Y, C);
    if true == PDE.show_subdomains
        PDE.showSubdomains('k', max(max(C)))
    end
    % zlim(Z_scale)
    view(0, 90)
    title(sprintf("Concentration\nNx=%d, Ny=%d, dt=%f, t=%f", Nx, Ny, dt, time))
    pause(0.00001)
end
hold off

end