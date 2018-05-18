function [ F ] = makeDirichletCondition( PDE, F, Nx, Ny, dx, dy )
% makeDirichletCondition: take into acount the Dirichlet condition
% specified in the PDE structure by modifying the F matrix
% [ F ] = makeDirichletCondition( PDE, F, Nx, Ny, dx, dy );
%
% Inputs:
%  PDE: PDE object containing information about the system
%  F:   pressure vector of size Ny by Nx, where Nx is the number of columns
%       and Ny the number of rows of the mesh
%  Nx:  number of discretization columns on the x axis
%  Ny:  number of discretization rows on the y axis
%  dx:  step size in the Ox direction
%  dy:  step size in the Oy direction
%
% Outputs:
%  F: pressure vector of size Ny by Nx, where Nx is the number of columns
%     and Ny the number of rows of the mesh


% Below and above borders
for i = 2 : Nx-1
    K = PDE.K.eval_matrix(PDE.x_min + dx*(i-1), PDE.y_min);
    gamma = 4*K(2,2)*dx^2;
    delta = (K(1,2) + K(2,1))*(dx*dy);
    F(1, i) = F(1, i) + delta * PDE.boundary_condition1(PDE.x_min + (i+1)*dx) + ...
                      - gamma * PDE.boundary_condition1(PDE.x_min + i*dx) + ...
                      - delta * PDE.boundary_condition1(PDE.x_min + (i-1)*dx);

    K = PDE.K.eval_matrix(PDE.x_min + dx*(i-1), PDE.y_max);
    gamma = 4*K(2,2)*dx^2;
    delta = (K(1,2) + K(2,1))*(dx*dy);
    F(Ny, i) = F(Ny, i) + delta * PDE.boundary_condition3(PDE.x_min + (i-1)*dx) + ...
                      - gamma * PDE.boundary_condition3(PDE.x_min + i*dx) + ...
                      - delta * PDE.boundary_condition3(PDE.x_min + (i+1)*dx);
end

% Left and right borders
for j = 2 : Ny-1
    K = PDE.K.eval_matrix(PDE.x_min, PDE.y_min + dy*(j-1));
    beta = 4*K(1,1)*dy^2;
    delta = (K(1,2) + K(2,1))*(dx*dy);
    F(j, 1) = F(j, 1) - delta * PDE.boundary_condition4(PDE.y_min + (j-1)*dy) + ...
                      - beta * PDE.boundary_condition4(PDE.y_min + j*dy) + ...
                      + delta * PDE.boundary_condition4(PDE.y_min + (j+1)*dy);

    beta = 4*PDE.K.eval_coef(1,1,PDE.x_max-dx,PDE.y_min + (j-1)*dy)*dy^2;
    delta =   (PDE.K.eval_coef(1,2,PDE.x_max-dx,PDE.y_min + (j-1)*dy) + ...
            + PDE.K.eval_coef(2,1,PDE.x_max-dx,PDE.y_min + (j-1)*dy))*(dx*dy);
    F(j, Nx) = F(j, Nx) - delta * PDE.boundary_condition2(PDE.y_min + (j+1)*dy) + ...
                        - beta * PDE.boundary_condition2(PDE.y_min + j*dy) + ...
                        + delta * PDE.boundary_condition2(PDE.y_min + (j-1)*dy);
end

% Lower left corner
K = PDE.K.eval_matrix(PDE.x_min, PDE.y_min);
beta = 4*K(1,1)*dy^2;
gamma = 4*K(2,2)*dx^2;
delta = (K(1,2) + K(2,1))*(dx*dy);
F(1,1) =   F(1,1) + ...
         + delta * PDE.boundary_condition1(PDE.x_min + 2*dx) ...
         - gamma * PDE.boundary_condition1(PDE.x_min + dx) ...
         - delta * PDE.boundary_condition4(PDE.y_min) ...
         - beta * PDE.boundary_condition4(PDE.y_min + dy) ...
         + delta * PDE.boundary_condition4(PDE.y_min + 2*dy);

% Lower right corner
K = PDE.K.eval_matrix(PDE.x_max, PDE.y_min);
beta = 4*K(1,1)*dy^2;
gamma = 4*K(2,2)*dx^2;
delta = (K(1,2) + K(2,1))*(dx*dy);
F(1, Nx) =   F(1, Nx) + ...
           - delta * PDE.boundary_condition2(PDE.y_min + 2*dy) ...
           - beta * PDE.boundary_condition2(PDE.y_min + dy) ...
           + delta * PDE.boundary_condition1(PDE.x_max) ...
           - gamma * PDE.boundary_condition1(PDE.x_max-dx) ...
           - delta * PDE.boundary_condition1(PDE.x_max-2*dx);

% Upper left corner
K = PDE.K.eval_matrix(PDE.x_min, PDE.y_max);
beta = 4*K(1,1)*dy^2;
gamma = 4*K(2,2)*dx^2;
delta = (K(1,2) + K(2,1))*(dx*dy);
F(Ny, 1) =   F(Ny, 1) + ...
           - delta * PDE.boundary_condition4(PDE.y_max-2*dy) ...
           - beta * PDE.boundary_condition4(PDE.y_max-dy) ...
           + delta * PDE.boundary_condition3(PDE.x_min) ...
           - gamma * PDE.boundary_condition3(PDE.x_min + dx) ...
           - delta * PDE.boundary_condition3(PDE.x_min + 2*dx);

% Upper right corner
K = PDE.K.eval_matrix(PDE.x_min, PDE.y_min);
beta = 4*K(1,1)*dy^2;
gamma = 4*K(2,2)*dx^2;
delta = (K(1,2) + K(2,1))*(dx*dy);
F(Ny,Nx) =   F(Ny,Nx) + ...
           + delta * PDE.boundary_condition3(PDE.x_max-2*dx) ...
           - gamma * PDE.boundary_condition3(PDE.x_max-dx) ...
           - delta * PDE.boundary_condition3(PDE.x_max) ...
           - beta * PDE.boundary_condition2(PDE.y_max-dy) ...
           + delta * PDE.boundary_condition2(PDE.y_max-2*dy);

end