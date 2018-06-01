clear
close all


%% PDE structure initialization

PDE = AdvectionDiffusionPDE();

PDE.x_min = -0.2; PDE.x_max = 1.0;
PDE.y_min = 0.0; PDE.y_max = 1.0;
PDE.T_final = .3;

PDE.Nx = 72;
PDE.Ny = 60;

[X, Y] = PDE.generateMesh();

% Second member (-div(K grad(h)) = f1)
f1_darcy1 =@(X, Y) zeros(size(X));    % Option 1
f1_darcy2 =@(X, Y) -4*ones(size(X));  % Option 2
f1_darcy3 =@(X, Y) -3*X - 6*Y.^2;     % Option 3
PDE.setSecondMemberDarcy(f1_darcy1);

% Initial condition
gaussian_curve =@(X, Y) .5*exp(-200 * ((X-.75).^2 + (Y-.75).^2));
PDE.setInitialCondition(gaussian_curve(X, Y));

% Dirichlet boundary conditions
exact_pressure1 =@(X, Y) 3+(2*X+2*Y);        % Option 1   
exact_pressure2 =@(X, Y) 3+(2*X.^2+2*Y.^2);  % Option 2
exact_pressure3 =@(X, Y) 3+(X.^2+Y.^2);      % Option 3
exact_pressure4 =@(X, Y) 1+0.25*(X.^3+Y.^4); % Option 4
PDE.setBoundaryConditions(exact_pressure1);

% Storing the permeability matrix
K = eye(2);
K(1,2) = .35;
K(2,1) = .35;
PDE.setPermeability(K);

% Storing the diffusion matrix
D = 0.02 * eye(2,2);
PDE.setDiffusion(D);
PDE.compute_diffusion = true;   % If disabled, no diffusion is taken into
                                % account and the program runs much faster


%% Geometry of the domain/obstacles

PDE.show_subdomains = false;    % If enabled, the domains are covered by a
                                % black surface

% First subdomain
A.x = 0.4; A.y = 0.4; B.x = 0.6; B.y = 0.6;  % Rectangle coordinates
subDomain1 = RectangleDomain(A, B, X, Y);    % Rectangle subdomain
K_subdomain1 = 0.00001 * eye(2,2);
D_subdomain1 = 0. * eye(2,2);
PDE.addSubdomain(subDomain1, K_subdomain1, D_subdomain1);

% Second subdomain
A.x = 0.3; A.y = 0.5; B.x = 0.35; B.y = 0.75;  % Rectangle coordinates
subDomain2 = RectangleDomain(A, B, X, Y);      % Rectangle subdomain
K_subdomain2 = 0.00001 * eye(2,2);
D_subdomain2 = 0. * eye(2,2);
PDE.addSubdomain(subDomain2, K_subdomain2, D_subdomain2);


%% Solving the system

H = solvePressure( PDE );
V = computeVelocity( PDE, H );

figure(1);
surf(X, Y, H)
hold on
surf(X, Y, reshape(V.x(:,2), [], PDE.Ny)')
title('Pressure and velocity (V.x, right borders)')
hold off

CFL = 1.0;

solveConcentration( PDE, V, CFL );
