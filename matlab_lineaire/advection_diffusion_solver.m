clear
close all


%% PDE structure initialization

PDE = AdvectionDiffusionPDE();

PDE.x_min = 0; PDE.x_max = 1.0;
PDE.y_min = 0; PDE.y_max = 1.0;
PDE.T_final = 0.5;

% Second members
f1_darcy1 =@(X, Y) -4*ones(size(X));
f1_darcy2 =@(X, Y) zeros(size(X));
PDE.setSecondMemberDarcy(f1_darcy1);

% Initial condition
gaussian_curve =@(X, Y) .5*exp(-100 * ((X-.75).^2 + (Y-.75).^2));
PDE.setInitialCondition(gaussian_curve);

% Dirichlet boundary conditions
exact_pressure1 =@(X, Y) 3+(2*X+2*Y);
exact_pressure2 =@(X, Y) 3+(2*X.^2+2*Y.^2);
exact_pressure3 =@(X, Y) 3+(X.^2+Y.^2);
PDE.setBoundaryConditions(exact_pressure3);

% Storing the permeability matrix
K = eye(2);
K(1,2) = .75;
PDE.setPermeability(K);

% Storing the diffusion matrix
D = zeros(2,2);
PDE.setDiffusion(D);


%% Geometry of the domain (obstacles)

PDE.show_subdomains = true;

% First subdomain
A.x = 0.4; A.y = 0.4; B.x = 0.6; B.y = 0.6;
subDomain1 = RectangleDomain(A, B);
K_subdomain1 = 0.001 * eye(2);
D_subdomain1 = D;
PDE.addSubdomain(subDomain1, K_subdomain1, D_subdomain1);


%% Domain discretization

Nx = 50;
Ny = 50;

[X, Y] = PDE.generateMesh(Nx, Ny);


%% Solving the system

h = solvePressure( PDE, Nx, Ny );
V = computeVelocity( PDE, h, Nx, Ny );

figure(1);
mesh(X, Y, h)
hold on
mesh(X, Y, reshape(V.x(:,2), [], Ny)')
title('Pressure and velocity (V.x, right borders)')
hold off

CFL = 0.5;
solveConcentration( PDE, V, Nx, Ny, CFL );


