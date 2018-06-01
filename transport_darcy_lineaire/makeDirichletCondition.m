function [ F ] = makeDirichletCondition( PDE, F, dx, dy )
% makeDirichletCondition: take into acount the Dirichlet condition
% specified in the PDE structure by modifying the F matrix
% [ F ] = makeDirichletCondition( PDE, F, Nx, Ny, dx, dy );
%
% Inputs:
%  PDE  PDE object containing information about the system
%  F -- pressure vector of size Ny by Nx, where Nx is the number of columns
%       and Ny the number of rows of the mesh
%  dx - step size in the Ox direction
%  dy - step size in the Oy direction
%
% Outputs:
%  F  pressure vector of size Ny by Nx, where Nx is the number of columns
%     and Ny the number of rows of the mesh

Nx = PDE.Nx-2;
Ny = PDE.Ny-2;

% Below and above borders
for i = 2 : Nx-1
    K12_left  = PDE.K.eval_harmonic_mean(1,2,i+1,2,[-1,0]);
    K12_right = PDE.K.eval_harmonic_mean(1,2,i+1,2,[1,0]);
    K21_below = PDE.K.eval_harmonic_mean(2,1,i+1,2,[0,-1]);
    K22_below = PDE.K.eval_harmonic_mean(2,2,i+1,2,[0,-1]);
    gamma = K12_left*dx*dy - K12_right*dx*dy + 4*K22_below*dx^2;
    delta1 = (K12_left + K21_below)*(dx*dy);
    delta2 = (K12_right + K21_below)*(dx*dy);
    F(1, i) = F(1, i) + delta2 * PDE.boundary_condition1(PDE.x_min + (i+1)*dx) + ...
                      - gamma * PDE.boundary_condition1(PDE.x_min + i*dx) + ...
                      - delta1 * PDE.boundary_condition1(PDE.x_min + (i-1)*dx);

    K12_left  = PDE.K.eval_harmonic_mean(1,2,i+1,Ny+1,[-1,0]);
    K12_right = PDE.K.eval_harmonic_mean(1,2,i+1,Ny+1,[1,0]);
    K21_above = PDE.K.eval_harmonic_mean(2,1,i+1,Ny+1,[0,1]);
    K22_above = PDE.K.eval_harmonic_mean(2,2,i+1,Ny+1,[0,1]);
    gamma = K12_right*dx*dy - K12_left*dx*dy + 4*K22_above*dx^2;
    delta1 = (K12_left + K21_above)*(dx*dy);
    delta2 = (K12_right + K21_above)*(dx*dy);
    F(Ny, i) = F(Ny, i) + delta1 * PDE.boundary_condition3(PDE.x_min + (i-1)*dx) + ...
                      - gamma * PDE.boundary_condition3(PDE.x_min + i*dx) + ...
                      - delta2 * PDE.boundary_condition3(PDE.x_min + (i+1)*dx);

end

% Left and right borders
for j = 2 : Ny-1
    K11_left  = PDE.K.eval_harmonic_mean(1,1,2,j+1,[-1,0]);
    K12_left  = PDE.K.eval_harmonic_mean(1,2,2,j+1,[-1,0]);
    K21_below = PDE.K.eval_harmonic_mean(2,1,2,j+1,[0,-1]);
    K21_above = PDE.K.eval_harmonic_mean(2,1,2,j+1,[0,1]);
    beta = K21_below*dx*dy - K21_above*dx*dy + 4*K11_left*dy^2;
    delta1 = (K12_left + K21_below)*(dx*dy);
    delta2 = (K12_left + K21_above)*(dx*dy);
    F(j, 1) = F(j, 1) - delta1 * PDE.boundary_condition4(PDE.y_min + (j-1)*dy) + ...
                      - beta * PDE.boundary_condition4(PDE.y_min + j*dy) + ...
                      + delta2 * PDE.boundary_condition4(PDE.y_min + (j+1)*dy);

    K11_right = PDE.K.eval_harmonic_mean(1,1,Nx+1,j+1,[1,0]);
    K12_right = PDE.K.eval_harmonic_mean(1,2,Nx+1,j+1,[1,0]);
    K21_below = PDE.K.eval_harmonic_mean(2,1,Nx+1,j+1,[0,-1]);
    K21_above = PDE.K.eval_harmonic_mean(2,1,Nx+1,j+1,[0,1]);
    beta = K21_above*dx*dy - K21_below*dx*dy + 4*K11_right*dy^2;
    delta1 = (K12_right + K21_above)*(dx*dy);
    delta2 = (K12_right + K21_below)*(dx*dy);
    F(j, Nx) = F(j, Nx) - delta1 * PDE.boundary_condition2(PDE.y_min + (j+1)*dy) + ...
                        - beta * PDE.boundary_condition2(PDE.y_min + j*dy) + ...
                        + delta2 * PDE.boundary_condition2(PDE.y_min + (j-1)*dy);
end

% Lower left corner
K11_left  = PDE.K.eval_harmonic_mean(1,1,2,2,[-1,0]);
K12_left  = PDE.K.eval_harmonic_mean(1,2,2,2,[-1,0]);
K12_right = PDE.K.eval_harmonic_mean(1,2,2,2,[1,0]);
K21_below = PDE.K.eval_harmonic_mean(2,1,2,2,[0,-1]);
K21_above = PDE.K.eval_harmonic_mean(2,1,2,2,[0,1]);
K22_below = PDE.K.eval_harmonic_mean(2,2,2,2,[0,-1]);
delta1 = (K12_right + K21_below)*(dx*dy);
gamma = K12_left*dx*dy - K12_right*dx*dy + 4*K22_below*dx^2;
delta2 = K12_left*dx*dy + K21_below*dx*dy;
beta = K21_below*dx*dy - K21_above*dx*dy + 4*K11_left*dy^2;
delta3 = K12_left*dx*dy + K21_above*dx*dy;
F(1,1) =   F(1,1) + ...
         + delta1 * PDE.boundary_condition1(PDE.x_min + 2*dx) ...
         - gamma * PDE.boundary_condition1(PDE.x_min + dx) ...
         - delta2 * PDE.boundary_condition4(PDE.y_min) ...
         - beta * PDE.boundary_condition4(PDE.y_min + dy) ...
         + delta3 * PDE.boundary_condition4(PDE.y_min + 2*dy);

% Lower right corner
K11_right = PDE.K.eval_harmonic_mean(1,1,Nx+1,2,[1,0]);
K12_left  = PDE.K.eval_harmonic_mean(1,2,Nx+1,2,[-1,0]);
K12_right = PDE.K.eval_harmonic_mean(1,2,Nx+1,2,[1,0]);
K21_below = PDE.K.eval_harmonic_mean(2,1,Nx+1,2,[0,-1]);
K21_above = PDE.K.eval_harmonic_mean(2,1,Nx+1,2,[0,1]);
K22_below = PDE.K.eval_harmonic_mean(2,2,Nx+1,2,[0,-1]);
delta1 = K12_left*dx*dy + K21_below*dx*dy;
gamma = K12_left*dx*dy - K12_right*dx*dy + 4*K22_below*dx^2;
delta2 = K12_right*dx*dy + K21_below*dx*dy;
beta = K21_above*dx*dy - K21_below*dx*dy + 4*K11_right*dy^2;
delta3 = K12_right*dx*dy + K21_above*dx*dy;
F(1, Nx) =   F(1, Nx) + ...
           - delta1 * PDE.boundary_condition1(PDE.x_max-2*dx) ...
           - gamma * PDE.boundary_condition1(PDE.x_max-dx) ...
           + delta2 * PDE.boundary_condition1(PDE.x_max) ...
           - beta * PDE.boundary_condition2(PDE.y_min + dy) ...
           - delta3 * PDE.boundary_condition2(PDE.y_min + 2*dy);

% Upper right corner
K11_right = PDE.K.eval_harmonic_mean(1,1,Nx+1,Ny+1,[1,0]);
K12_left  = PDE.K.eval_harmonic_mean(1,2,Nx+1,Ny+1,[-1,0]);
K12_right = PDE.K.eval_harmonic_mean(1,2,Nx+1,Ny+1,[1,0]);
K21_below = PDE.K.eval_harmonic_mean(2,1,Nx+1,Ny+1,[0,-1]);
K21_above = PDE.K.eval_harmonic_mean(2,1,Nx+1,Ny+1,[0,1]);
K22_above = PDE.K.eval_harmonic_mean(2,2,Nx+1,Ny+1,[0,1]);
delta1 = K12_right*dx*dy + K21_below*dx*dy;
beta = K21_above*dx*dy - K21_below*dx*dy + 4*K11_right*dy^2;
delta2 = K12_right*dx*dy + K21_above*dx*dy;
gamma = K12_right*dx*dy - K12_left*dx*dy + 4*K22_above*dx^2;
delta3 = K12_left*dx*dy + K21_above*dx*dy;
F(Ny,Nx) =   F(Ny,Nx) + ...
           + delta1 * PDE.boundary_condition2(PDE.y_max-2*dy) ...
           - beta * PDE.boundary_condition2(PDE.y_max-dy) ...
           - delta2 * PDE.boundary_condition3(PDE.x_max) ...
           - gamma * PDE.boundary_condition3(PDE.x_max-dx) ...
           + delta3 * PDE.boundary_condition3(PDE.x_max-2*dx);

% Upper left corner
K11_left  = PDE.K.eval_harmonic_mean(1,1,2,Ny+1,[-1,0]);
K12_left  = PDE.K.eval_harmonic_mean(1,2,2,Ny+1,[-1,0]);
K12_right = PDE.K.eval_harmonic_mean(1,2,2,Ny+1,[1,0]);
K21_below = PDE.K.eval_harmonic_mean(2,1,2,Ny+1,[0,-1]);
K21_above = PDE.K.eval_harmonic_mean(2,1,2,Ny+1,[0,1]);
K22_above = PDE.K.eval_harmonic_mean(2,2,2,Ny+1,[0,1]);
delta1 = K12_right*dx*dy + K21_above*dx*dy;
gamma = K12_right*dx*dy - K12_left*dx*dy + 4*K22_above*dx^2;
delta2 = K12_left*dx*dy + K21_above*dx*dy;
beta = K21_below*dx*dy - K21_above*dx*dy + 4*K11_left*dy^2;
delta3 = K12_left*dx*dy + K21_below*dx*dy;
F(Ny, 1) =   F(Ny, 1) + ...
           - delta1 * PDE.boundary_condition3(PDE.x_min + 2*dx) ...
           - gamma * PDE.boundary_condition3(PDE.x_min + dx) ...
           + delta2 * PDE.boundary_condition3(PDE.x_min) ...
           - beta * PDE.boundary_condition4(PDE.y_max-dy) ...
           - delta3 * PDE.boundary_condition4(PDE.y_max-2*dy);

end