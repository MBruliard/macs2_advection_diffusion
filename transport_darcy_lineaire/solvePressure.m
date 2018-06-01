function [ H ] = solvePressure( PDE )
% solvePressure: return the pressure matrix H given the PDE system.
% [ H ] = solvePressure( PDE, Nx, Ny );
%
% Inputs:
%  PDE  PDE object containing information about the system
%
% Outputs:
%  H  pressure matrix of size Ny by Nx, where Nx is the number of columns
%     and Ny the number of rows of the mesh


dx = (PDE.x_max - PDE.x_min) / (PDE.Nx-1);
dy = (PDE.y_max - PDE.y_min) / (PDE.Ny-1);

% Rigidity matrix
R = vectorizedRigidityMatrix( PDE, dx, dy );

% Second member definition
[X, Y] = meshgrid(PDE.x_min + dx : dx : PDE.x_max - dx, ...
    PDE.y_min + dy : dy : PDE.y_max - dy);
F_matrix = -4*dx^2*dy^2 * PDE.f1(X, Y);
F_matrix = makeDirichletCondition( PDE, F_matrix, dx, dy );

% We reshape the matrix F_matrix into a vector
F_vector = reshape(F_matrix', [], 1);

% We compute the pressure
h = R \ F_vector;

% We reshpe the vector h into matrix
H = zeros(PDE.Ny, PDE.Nx);

% We set the values of H in the inner part of the mesh
H(2:PDE.Ny-1, 2:PDE.Nx-1) = reshape(h, [], PDE.Ny-2)';

% We set the values of H on the lower and upper borders
y_coord = 1 : PDE.Ny;
H(y_coord, 1) = PDE.boundary_condition4(PDE.y_min + (y_coord-1)*dy);
H(y_coord, PDE.Nx) = PDE.boundary_condition2(PDE.y_min + (y_coord-1)*dy);

% We set the values of H on the left and right borders
x_coord = 1 : PDE.Nx;
H(1, x_coord) = PDE.boundary_condition1(PDE.x_min + (x_coord-1)*dx);
H(PDE.Ny, x_coord) = PDE.boundary_condition3(PDE.x_min + (x_coord-1)*dx);

end