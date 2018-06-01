function [ C ] = solveConcentration( PDE, V, CFL )
% solveConcentration: approximates the concentration matrix C given the PDE
% system between times t=0 and t=PDE.T_final, and returns its last state.
% [ C ] = solveConcentration( PDE, V, Nx, Ny );
%
% Inputs:
%  PDE  PDE object containing information about the system
%  V -- velocity object with two subfields V.x and V.y
%  V.x  velocity matrix (Ox direction) of size Nx x Ny by 4, where the
%       first argument is the cell number and the second one is the
%       border index (anticlockwise, lower border is number 1)
%  V.y  velocity matrix (Oy direction) of size Nx x Ny by 4, where the
%       first argument is the cell number and the second one is the
%       border index (anticlockwise, lower border is number 1)
%  CFL  float value bounded by 0 and 1


dx = (PDE.x_max - PDE.x_min) / (PDE.Nx-1);
dy = (PDE.y_max - PDE.y_min) / (PDE.Ny-1);

nb_subdomains = length(PDE.subdomains);

% Initial conditions
[X, Y] = PDE.generateMesh();
C = PDE.initial_condition;

% Time step (CFL consition)
Vx_max = max(max(abs(V.x)));
Vy_max = max(max(abs(V.y)));
dt_advection = dx*dy / (dy*Vx_max + dx*Vy_max);

m = PDE.D.get_max_coef(1:PDE.Nx, 1:PDE.Ny);
if true == PDE.compute_diffusion && 0 ~= m
    dt_diffusion = 0.5 * dx^2*dy^2 / ((dx^2+dy^2) * m);
else
    dt_diffusion = dt_advection;
    m = 0;
end

dt = CFL * 0.5 * min(dt_advection, dt_diffusion);

title_format_string = 'Concentration evolution (linear solver), C_{total}=%f\nNx=%d, Ny=%d, dt=%f, time=%f, Nb subdomains=%d\nmax(Vx)=%f, max(Vy)=%f, max diffusion=%f';

hold on
% We loop until final time is reached
time = 0;
while time < PDE.T_final
    dt_ = min(dt, PDE.T_final-time);
    time = time + dt_;
    C = updateConcentration( PDE, C, V, dt_ );
    C_tot = sum(sum(C));
    figure(2);
    surf(X, Y, C, 'edgecolor', 'none');
    if true == PDE.show_subdomains
        PDE.showSubdomains('k', C_max)
    end
    view(0, 90)
    title(sprintf(title_format_string, C_tot, PDE.Nx, PDE.Ny, dt, time, nb_subdomains, Vx_max, Vy_max, m))
    set(gcf, 'Position', [500, 300, 650, 500], 'units','pix')
    set(gca, 'TickDir', 'out', 'TickLength', [.02 .02], ...
             'XMinorTick', 'on', 'YMinorTick', 'on', 'LineWidth', 1)
    pause(0.)

end
hold off


end