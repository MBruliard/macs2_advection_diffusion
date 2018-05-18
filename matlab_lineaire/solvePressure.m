function [ H ] = solvePressure( PDE, Nx, Ny )
% solvePressure: return the pressure matrix H given the PDE system.
% [ H ] = solvePressure( PDE, Nx, Ny );
%
% Inputs:
%  PDE  PDE object containing information about the system
%  Nx - number of discretization columns on the x axis
%  Ny - number of discretization rows on the y axis
%
% Outputs:
%  H  pressure matrix of size Ny by Nx, where Nx is the number of columns
%     and Ny the number of rows of the mesh


dx = (PDE.x_max - PDE.x_min) / (Nx-1);
dy = (PDE.x_max - PDE.x_min) / (Ny-1);

% Rigidity matrix
R = rigidityMatrix( PDE, Nx-2, Ny-2, dx, dy );

% Second member definition
[X, Y] = meshgrid(PDE.x_min + dx : dx : PDE.x_max - dx, ...
                  PDE.y_min + dy : dy : PDE.y_max - dy);
F_matrix = -8*dx^2*dy^2 * PDE.f1(X, Y);
F_matrix = makeDirichletCondition( PDE, F_matrix, Nx-2, Ny-2, dx, dy );

% We reshape the matrix F_matrix into a vector
F_vector = zeros((Ny-2)*(Nx-2), 1);
for j = 1 : Ny-2
    for i = 1 : Nx-2
        k = i + (j-1)*(Nx-2);
        F_vector(k) = F_matrix(j, i);
    end
end

% We compute the pressure
h = R \ F_vector;

% We reshpe the vector h into matrix
H = zeros(Ny, Nx);

% We set the values of H in the inner part of the mesh
for j = 2 : Ny-1
    for i = 2 : Nx-1
        k = i-1 + (j-2)*(Nx-2);
        H(j, i) = h(k);
    end
end

% We set the values of H on the borders
% Lower and upper borders
for i = 1 : Ny
    H(i, 1) = PDE.boundary_condition4(PDE.y_min + (i-1)*dy);
    H(i, Nx) = PDE.boundary_condition2(PDE.y_min + (i-1)*dy);
end
% Left and right borders
for j = 1 : Nx
    H(1, j) = PDE.boundary_condition1(PDE.x_min + (j-1)*dx);
    H(Ny, j) = PDE.boundary_condition3(PDE.x_min + (j-1)*dx);
end

end