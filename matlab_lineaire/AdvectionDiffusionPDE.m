classdef AdvectionDiffusionPDE < handle
    % AdvectionDiffusionPDE: class containing information about the PDE
    % system
    
    %% PROPERTIES DECLARATION

    properties % (Access = public)
        x_min
        x_max
        y_min
        y_max
        T_final
        K
        D
        f1
        f2
        initial_condition
        boundary_condition1
        boundary_condition2
        boundary_condition3
        boundary_condition4
        subdomains
        show_subdomains
        compute_diffusion
    end
    
    
    %% METHODS DEFINITION

    methods
        % Constructor
        function obj = AdvectionDiffusionPDE()
            obj.show_subdomains = true;
        end
        
        % Storing the permeability matrix
        function obj = setPermeability(obj, K)
            obj.K = MatrixFunction(K);
        end
        
        % Storing the diffusion matrix
        function obj = setDiffusion(obj, D)
            obj.D = MatrixFunction(D);
        end
        
        % Setting the boundary conditions as the trace of the given
        % function
        function obj = setBoundaryConditions(obj, f)
            obj.boundary_condition1 =@(X) f(X, obj.y_min);
            obj.boundary_condition2 =@(Y) f(obj.x_max, Y);
            obj.boundary_condition3 =@(X) f(X, obj.y_max);
            obj.boundary_condition4 =@(Y) f(obj.x_min, Y);
        end
        
        % Setting the second member occuring in the pressure equation
        function obj = setSecondMemberDarcy(obj, f)
            obj.f1 = f;
        end
        
        % Setting the second member occuring in the convection equation
        function obj = setSecondMemberConvection(obj, f)
            obj.f2 = f;
        end
        
        % Setting the initial condition
        function obj = setInitialCondition(obj, f)
            obj.initial_condition = f;
        end
        
        % Add a subdomain by giving its gemeometry and the associated
        % values of the permeability and diffusion matrix
        function obj = addSubdomain(obj, domain, K, D)
            nb_subdomains = length(obj.subdomains);
            obj.subdomains{nb_subdomains+1} = domain;
            obj.K.addSubdomain(domain, K);
            obj.D.addSubdomain(domain, D);
        end
        
        % Plot the subdomain geometry at position z_pos on the Z-axis
        function showSubdomains(obj, color, z_pos)
            nb_subdomains = length(obj.subdomains);
            for i = 1 : nb_subdomains
                obj.subdomains{i}.show(color, z_pos);
            end
        end
        
        % Generate a mesh
        function [X, Y] = generateMesh(obj, Nx, Ny)
            dx = (obj.x_max - obj.x_min) / (Nx-1);
            dy = (obj.y_max - obj.y_min) / (Ny-1);
            [X, Y] = meshgrid(obj.x_min : dx : obj.x_max, ...
                              obj.y_min : dy : obj.y_max);
        end
    end

end

