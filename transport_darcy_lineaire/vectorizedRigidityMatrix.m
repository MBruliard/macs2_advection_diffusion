function [ R ] = vectorizedRigidityMatrix( PDE, dx, dy )
% vectorizedRigidityMatrix: return the rigidity matrix R.
% [ R ] = vectorizedRigidityMatrix( PDE, dx, dy );
%
% Inputs:
%  PDE  PDE object containing information about the system
%  dx - step size in the Ox direction
%  dy - step size in the Oy direction
%
% Outputs:
%  R  matrix of size Ny x Nx by Ny x Nx, where Nx is the number of columns
%     and Ny the number of rows of the mesh


Nx = PDE.Nx-2;
Ny = PDE.Ny-2;

nb_nonzero = (3*(Ny-2) + 4) * (3*(Nx-2) + 4);

row_indexes    = zeros(1,nb_nonzero);
column_indexes = zeros(1,nb_nonzero);
nonzero_values = zeros(1,nb_nonzero);
vect_pos       = 1;


%% INTERIOR MESH

nx_ = 2 : Nx-1;
x_index_ = nx_+1;

for ny = 2 : Ny-1

    y_index = ny+1;
    K_12_left_  = PDE.K.eval_harmonic_mean_(1,2,x_index_,y_index,[-1, 0]);
    K_21_below_ = PDE.K.eval_harmonic_mean_(2,1,x_index_,y_index,[0, -1]);
    K_12_right_ = PDE.K.eval_harmonic_mean_(1,2,x_index_,y_index,[1, 0]);
    K_22_below_ = PDE.K.eval_harmonic_mean_(2,2,x_index_,y_index,[0, -1]);
    K_21_above_ = PDE.K.eval_harmonic_mean_(2,1,x_index_,y_index,[0, 1]);
    K_11_left_  = PDE.K.eval_harmonic_mean_(1,1,x_index_,y_index,[-1, 0]);
    K_11_right_ = PDE.K.eval_harmonic_mean_(1,1,x_index_,y_index,[1, 0]);
    K_22_above_ = PDE.K.eval_harmonic_mean_(2,2,x_index_,y_index,[0, 1]);
    
    k_ = nx_ + (ny-1) * Nx;
    
    % Coefficient (k,k-Nx-1)
    delta_ = K_12_left_*dx*dy + K_21_below_*dx*dy;
    row_indexes(vect_pos:vect_pos+Nx-3) = k_;
    column_indexes(vect_pos:vect_pos+Nx-3) = k_-Nx-1;
    nonzero_values(vect_pos:vect_pos+Nx-3) = delta_;
    vect_pos = vect_pos+Nx-2;
    
    % Coefficient (k,k-Nx)
    gamma_ = K_12_left_*dx*dy - K_12_right_*dx*dy + 4*K_22_below_*dx^2;
    row_indexes(vect_pos:vect_pos+Nx-3) = k_;
    column_indexes(vect_pos:vect_pos+Nx-3) = k_-Nx;
    nonzero_values(vect_pos:vect_pos+Nx-3) = gamma_;
    vect_pos = vect_pos+Nx-2;
    
    % Coefficient (k,k-Nx+1)
    delta_ = K_12_right_*dx*dy + K_21_below_*dx*dy;
    row_indexes(vect_pos:vect_pos+Nx-3) = k_;
    column_indexes(vect_pos:vect_pos+Nx-3) = k_-Nx+1;
    nonzero_values(vect_pos:vect_pos+Nx-3) = -delta_;
    vect_pos = vect_pos+Nx-2;
    
    % Coefficient (k,k-1)
    beta_ = K_21_below_*dx*dy - K_21_above_*dx*dy + 4*K_11_left_*dy^2;
    row_indexes(vect_pos:vect_pos+Nx-3) = k_;
    column_indexes(vect_pos:vect_pos+Nx-3) = k_-1;
    nonzero_values(vect_pos:vect_pos+Nx-3) = beta_;
    vect_pos = vect_pos+Nx-2;
    
    % Coefficient (k,k)
    alpha_ = - 4*K_11_left_*dy^2  - 4*K_11_right_*dy^2 + ...
             - 4*K_22_below_*dx^2 - 4*K_22_above_*dx^2;
    row_indexes(vect_pos:vect_pos+Nx-3) = k_;
    column_indexes(vect_pos:vect_pos+Nx-3) = k_;
    nonzero_values(vect_pos:vect_pos+Nx-3) = alpha_;
    vect_pos = vect_pos+Nx-2;
    
    % Coefficient (k,k+1)
    beta_ = K_21_above_*dx*dy - K_21_below_*dx*dy + 4*K_11_right_*dy^2;
    row_indexes(vect_pos:vect_pos+Nx-3) = k_;
    column_indexes(vect_pos:vect_pos+Nx-3) = k_+1;
    nonzero_values(vect_pos:vect_pos+Nx-3) = beta_;
    vect_pos = vect_pos+Nx-2;
    
    % Coefficient (k,k+Nx-1)
    delta_ = K_12_left_*dx*dy + K_21_above_*dx*dy;
    row_indexes(vect_pos:vect_pos+Nx-3) = k_;
    column_indexes(vect_pos:vect_pos+Nx-3) = k_+Nx-1;
    nonzero_values(vect_pos:vect_pos+Nx-3) = -delta_;
    vect_pos = vect_pos+Nx-2;

    % Coefficient (k,k+Nx)
    gamma_ = K_12_right_*dx*dy - K_12_left_*dx*dy + 4*K_22_above_*dx^2;
    row_indexes(vect_pos:vect_pos+Nx-3) = k_;
    column_indexes(vect_pos:vect_pos+Nx-3) = k_+Nx;
    nonzero_values(vect_pos:vect_pos+Nx-3) = gamma_;
    vect_pos = vect_pos+Nx-2;

    % Coefficient (k,k+Nx+1)
    delta_ = K_12_right_*dx*dy + K_21_above_*dx*dy;
    row_indexes(vect_pos:vect_pos+Nx-3) = k_;
    column_indexes(vect_pos:vect_pos+Nx-3) = k_+Nx+1;
    nonzero_values(vect_pos:vect_pos+Nx-3) = delta_;
    vect_pos = vect_pos+Nx-2;
    
end


%% LOWER BORDER MESH

ny = 1;
for nx = 2 : Nx-1
    
    k = nx + (ny-1) * Nx;
    x_index = nx+1;
    y_index = ny+1;
    
    K_21_below = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, -1]);
    K_21_above = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, 1]);
    K_11_left  = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[-1, 0]);
    K_11_right = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[1, 0]);
    K_22_below = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, -1]);
    K_22_above = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, 1]);
    K_12_left  = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[-1, 0]);
    K_12_right = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[1, 0]);
    
    % Coefficient (k,k-1)
    beta = K_21_below*dx*dy - K_21_above*dx*dy + 4*K_11_left*dy^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k-1;
    nonzero_values(vect_pos) = beta;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k)
    alpha = - 4*K_11_left*dy^2  - 4*K_11_right*dy^2 + ...
            - 4*K_22_below*dx^2 - 4*K_22_above*dx^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k;
    nonzero_values(vect_pos) = alpha;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k+1)
    beta = K_21_above*dx*dy - K_21_below*dx*dy + 4*K_11_right*dy^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k+1;
    nonzero_values(vect_pos) = beta;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k+Nx-1)
    delta = K_12_left*dx*dy + K_21_above*dx*dy;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k+Nx-1;
    nonzero_values(vect_pos) = -delta;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k+Nx)
    gamma = K_12_right*dx*dy - K_12_left*dx*dy + 4*K_22_above*dx^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k+Nx;
    nonzero_values(vect_pos) = gamma;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k+Nx+1)
    delta = K_12_right*dx*dy + K_21_above*dx*dy;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k+Nx+1;
    nonzero_values(vect_pos) = delta;
    vect_pos = vect_pos+1;
    
end


%% UPPER BORDER MESH

ny = Ny;
for nx = 2 : Nx-1
    
    k = nx + (ny-1) * Nx;
    x_index = nx+1;
    y_index = ny+1;
    
    K_12_left  = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[-1, 0]);
    K_12_right = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[1, 0]);
    K_11_left  = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[-1, 0]);
    K_22_below = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, -1]);
    K_22_above = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, 1]);
    K_21_above = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, 1]);
    K_21_below = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, -1]);
    K_11_right = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[1, 0]);

    % Coefficient (k,k-Nx-1)
    delta = K_12_left*dx*dy + K_21_below*dx*dy;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k-Nx-1;
    nonzero_values(vect_pos) = delta;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k-Nx)
    gamma = K_12_left*dx*dy - K_12_right*dx*dy + 4*K_22_below*dx^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k-Nx;
    nonzero_values(vect_pos) = gamma;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k-Nx+1)
    delta = K_12_right*dx*dy + K_21_below*dx*dy;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k-Nx+1;
    nonzero_values(vect_pos) = -delta;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k-1)
    beta = K_21_below*dx*dy - K_21_above*dx*dy + 4*K_11_left*dy^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k-1;
    nonzero_values(vect_pos) = beta;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k)
    alpha = - 4*K_11_left*dy^2  - 4*K_11_right*dy^2 + ...
            - 4*K_22_below*dx^2 - 4*K_22_above*dx^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k;
    nonzero_values(vect_pos) = alpha;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k+1)
    beta = K_21_above*dx*dy - K_21_below*dx*dy + 4*K_11_right*dy^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k+1;
    nonzero_values(vect_pos) = beta;
    vect_pos = vect_pos+1;
    
end


%% LEFT BORDER MESH

nx = 1;
for ny = 2 : Ny-1

    k = nx + (ny-1) * Nx;
    x_index = nx+1;
    y_index = ny+1;

    K_11_left  = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[-1, 0]);
    K_22_below = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, -1]);
    K_21_below = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, -1]);
    K_11_right = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[1, 0]);
    K_12_left  = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[-1, 0]);
    K_22_above = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, 1]);
    K_12_right = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[1, 0]);
    K_21_above = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, 1]);
    
    % Coefficient (k,k-Nx)
    gamma = K_12_left*dx*dy - K_12_right*dx*dy + 4*K_22_below*dx^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k-Nx;
    nonzero_values(vect_pos) = gamma;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k-Nx+1)
    delta = K_12_right*dx*dy + K_21_below*dx*dy;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k-Nx+1;
    nonzero_values(vect_pos) = -delta;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k)
    alpha = - 4*K_11_left*dy^2  - 4*K_11_right*dy^2 + ...
            - 4*K_22_below*dx^2 - 4*K_22_above*dx^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k;
    nonzero_values(vect_pos) = alpha;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k+1)
    beta = K_21_above*dx*dy - K_21_below*dx*dy + 4*K_11_right*dy^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k+1;
    nonzero_values(vect_pos) = beta;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k+Nx)
    gamma = K_12_right*dx*dy - K_12_left*dx*dy + 4*K_22_above*dx^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k+Nx;
    nonzero_values(vect_pos) = gamma;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k+Nx+1)
    delta = K_12_right*dx*dy + K_21_above*dx*dy;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k+Nx+1;
    nonzero_values(vect_pos) = delta;
    vect_pos = vect_pos+1;
    
end


%% RIGHT BORDER MESH

nx = Nx;
for ny = 2 : Ny-1
    
    k = nx + (ny-1) * Nx;
    x_index = nx+1;
    y_index = ny+1;
    
    K_21_below = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, -1]);
    K_11_left  = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[-1, 0]);
    K_11_right = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[1, 0]);
    K_22_below = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, -1]);
    K_21_above = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, 1]);
    K_12_right = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[1, 0]);
    K_12_left  = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[-1, 0]);
    K_22_above = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, 1]);

    % Coefficient (k,k-Nx-1)
    delta = K_12_left*dx*dy + K_21_below*dx*dy;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k-Nx-1;
    nonzero_values(vect_pos) = delta;
    vect_pos = vect_pos+1;

    % Coefficient (k,k-Nx)
    gamma = K_12_left*dx*dy - K_12_right*dx*dy + 4*K_22_below*dx^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k-Nx;
    nonzero_values(vect_pos) = gamma;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k-1)
    beta = K_21_below*dx*dy - K_21_above*dx*dy + 4*K_11_left*dy^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k-1;
    nonzero_values(vect_pos) = beta;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k)
    alpha = - 4*K_11_left*dy^2  - 4*K_11_right*dy^2 + ...
            - 4*K_22_below*dx^2 - 4*K_22_above*dx^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k;
    nonzero_values(vect_pos) = alpha;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k+Nx-1)
    delta = K_12_left*dx*dy + K_21_above*dx*dy;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k+Nx-1;
    nonzero_values(vect_pos) = -delta;
    vect_pos = vect_pos+1;
    
    % Coefficient (k,k+Nx)
    gamma = K_12_right*dx*dy - K_12_left*dx*dy + 4*K_22_above*dx^2;
    row_indexes(vect_pos) = k;
    column_indexes(vect_pos) = k+Nx;
    nonzero_values(vect_pos) = gamma;
    vect_pos = vect_pos+1;
    
end


%% LOWER LEFT CELL

nx = 1;
ny = 1;
k = nx + (ny-1) * Nx;
x_index = nx+1;
y_index = ny+1;


K_11_left  = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[-1, 0]);
K_22_below = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, -1]);
K_21_below = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, -1]);
K_11_right = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[1, 0]);
K_12_left  = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[-1, 0]);
K_22_above = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, 1]);
K_12_right = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[1, 0]);
K_21_above = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, 1]);

% Coefficient (k,k)
alpha = - 4*K_11_left*dy^2  - 4*K_11_right*dy^2 + ...
        - 4*K_22_below*dx^2 - 4*K_22_above*dx^2;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k;
nonzero_values(vect_pos) = alpha;
vect_pos = vect_pos+1;

% Coefficient (k,k+1)
beta = K_21_above*dx*dy - K_21_below*dx*dy + 4*K_11_right*dy^2;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k+1;
nonzero_values(vect_pos) = beta;
vect_pos = vect_pos+1;

% Coefficient (k,k+Nx)
gamma = K_12_right*dx*dy - K_12_left*dx*dy + 4*K_22_above*dx^2;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k+Nx;
nonzero_values(vect_pos) = gamma;
vect_pos = vect_pos+1;

% Coefficient (k,k+Nx+1)
delta = K_12_right*dx*dy + K_21_above*dx*dy;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k+Nx+1;
nonzero_values(vect_pos) = delta;
vect_pos = vect_pos+1;


%% LOWER RIGHT CELL

nx = Nx;
ny = 1;
k = nx + (ny-1) * Nx;
x_index = nx+1;
y_index = ny+1;

K_21_below = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, -1]);
K_11_left  = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[-1, 0]);
K_11_right = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[1, 0]);
K_22_below = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, -1]);
K_21_above = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, 1]);
K_12_right = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[1, 0]);
K_12_left  = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[-1, 0]);
K_22_above = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, 1]);

% Coefficient (k,k-1)
beta = K_21_below*dx*dy - K_21_above*dx*dy + 4*K_11_left*dy^2;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k-1;
nonzero_values(vect_pos) = beta;
vect_pos = vect_pos+1;

% Coefficient (k,k)
alpha = - 4*K_11_left*dy^2  - 4*K_11_right*dy^2 + ...
        - 4*K_22_below*dx^2 - 4*K_22_above*dx^2;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k;
nonzero_values(vect_pos) = alpha;
vect_pos = vect_pos+1;

% Coefficient (k,k+Nx-1)
delta = K_12_left*dx*dy + K_21_above*dx*dy;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k+Nx-1;
nonzero_values(vect_pos) = -delta;
vect_pos = vect_pos+1;

% Coefficient (k,k+Nx)
gamma = K_12_right*dx*dy - K_12_left*dx*dy + 4*K_22_above*dx^2;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k+Nx;
nonzero_values(vect_pos) = gamma;
vect_pos = vect_pos+1;


%% UPPER RIGHT CELL

nx = Nx;
ny = Ny;
k = nx + (ny-1) * Nx;
x_index = nx+1;
y_index = ny+1;

K_12_left  = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[-1, 0]);
K_12_right = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[1, 0]);
K_21_below = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, -1]);
K_21_above = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, 1]);
K_11_left  = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[-1, 0]);
K_11_right = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[1, 0]);
K_22_below = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, -1]);
K_22_above = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, 1]);

% Coefficient (k,k-Nx-1)
delta = K_12_left*dx*dy + K_21_below*dx*dy;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k-Nx-1;
nonzero_values(vect_pos) = delta;
vect_pos = vect_pos+1;

% Coefficient (k,k-Nx)
gamma = K_12_left*dx*dy - K_12_right*dx*dy + 4*K_22_below*dx^2;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k-Nx;
nonzero_values(vect_pos) = gamma;
vect_pos = vect_pos+1;

% Coefficient (k,k-1)
beta = K_21_below*dx*dy - K_21_above*dx*dy + 4*K_11_left*dy^2;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k-1;
nonzero_values(vect_pos) = beta;
vect_pos = vect_pos+1;

% Coefficient (k,k)
alpha = - 4*K_11_left*dy^2  - 4*K_11_right*dy^2 + ...
        - 4*K_22_below*dx^2 - 4*K_22_above*dx^2;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k;
nonzero_values(vect_pos) = alpha;
vect_pos = vect_pos+1;


%% UPPER LEFT CELL

nx = 1;
ny = Ny;
k = nx + (ny-1) * Nx;
x_index = nx+1;
y_index = ny+1;

K_12_left  = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[-1, 0]);
K_12_right = PDE.K.eval_harmonic_mean(1,2,x_index,y_index,[1, 0]);
K_11_left  = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[-1, 0]);
K_22_below = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, -1]);
K_22_above = PDE.K.eval_harmonic_mean(2,2,x_index,y_index,[0, 1]);
K_21_above = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, 1]);
K_21_below = PDE.K.eval_harmonic_mean(2,1,x_index,y_index,[0, -1]);
K_11_right = PDE.K.eval_harmonic_mean(1,1,x_index,y_index,[1, 0]);

% Coefficient (k,k-Nx)
gamma = K_12_left*dx*dy - K_12_right*dx*dy + 4*K_22_below*dx^2;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k-Nx;
nonzero_values(vect_pos) = gamma;
vect_pos = vect_pos+1;

% Coefficient (k,k-Nx+1)
delta = K_12_right*dx*dy + K_21_below*dx*dy;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k-Nx+1;
nonzero_values(vect_pos) = -delta;
vect_pos = vect_pos+1;

% Coefficient (k,k)
alpha = - 4*K_11_left*dy^2  - 4*K_11_right*dy^2 + ...
        - 4*K_22_below*dx^2 - 4*K_22_above*dx^2;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k;
nonzero_values(vect_pos) = alpha;
vect_pos = vect_pos+1;

% Coefficient (k,k+1)
beta = K_21_above*dx*dy - K_21_below*dx*dy + 4*K_11_right*dy^2;
row_indexes(vect_pos) = k;
column_indexes(vect_pos) = k+1;
nonzero_values(vect_pos) = beta;


%% ASSEMBLING SPARSE RIGIDITY MATRIX

R = sparse(row_indexes, column_indexes, nonzero_values, Ny*Nx, Ny*Nx);

end