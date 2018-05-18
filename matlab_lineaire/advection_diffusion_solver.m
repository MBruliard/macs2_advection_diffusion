clear
close all


%% PDE structure initialization

PDE = AdvectionDiffusionPDE();

PDE.x_min = 0; PDE.x_max = 1.0;
PDE.y_min = 0; PDE.y_max = 1.0;
PDE.T_final = 0.4;

% Second member (-div(K grad(h)) = f1)
f1_darcy1 =@(X, Y) zeros(size(X));    % Option 1
f1_darcy2 =@(X, Y) -4*ones(size(X));  % Option 2
PDE.setSecondMemberDarcy(f1_darcy1);

% Initial condition
gaussian_curve =@(X, Y) .5*exp(-100 * ((X-.75).^2 + (Y-.75).^2));
PDE.setInitialCondition(gaussian_curve);

% Dirichlet boundary conditions
exact_pressure1 =@(X, Y) 3+(2*X+2*Y);        % Option 1   
exact_pressure2 =@(X, Y) 3+(2*X.^2+2*Y.^2);  % Option 2
exact_pressure3 =@(X, Y) 3+(X.^2+Y.^2);      % Option 3
PDE.setBoundaryConditions(exact_pressure1);

% Storing the permeability matrix
K = 1 * eye(2);
% K(1,2) = .75;
PDE.setPermeability(K);

% Storing the diffusion matrix
D = 0.05 * eye(2,2);
PDE.setDiffusion(D);
PDE.compute_diffusion = true;   % If disabled, no diffusion is taken into
                                % account and the program runs much faster


%% Geometry of the domain/obstacles

PDE.show_subdomains = false;    % If enabled, the domains are covered by a
                                % black surface

% First subdomain
A.x = 0.4; A.y = 0.4; B.x = 0.6; B.y = 0.6;  % Rectangle coordinates
subDomain1 = RectangleDomain(A, B);          % Rectangle subdomain
K_subdomain1 = 0.0001 * eye(2,2);
D_subdomain1 = 0.0000 * eye(2,2);
PDE.addSubdomain(subDomain1, K_subdomain1, D_subdomain1);


%% Domain discretization

Nx = 52;
Ny = 52;

[X, Y] = PDE.generateMesh(Nx, Ny);


%% Solving the system

H = solvePressure( PDE, Nx, Ny );
V = computeVelocity( PDE, H, Nx, Ny );

figure(1);
mesh(X, Y, H)
hold on
mesh(X, Y, reshape(V.x(:,2), [], Ny)')
title('Pressure and velocity (V.x, right borders)')
hold off

CFL = 1.0;
solveConcentration( PDE, V, Nx, Ny, CFL );


