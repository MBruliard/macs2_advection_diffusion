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
        Nx
        Ny
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
        % CONSTRUCTOR
        function obj = AdvectionDiffusionPDE()
            obj.show_subdomains = true;
        end
        
        % STORING THE PERMEABILITY MATRIX
        function obj = setPermeability(obj, K)
            obj.K = MatrixFunction(K, obj.Nx, obj.Ny);
        end
        
        % STORING THE DIFFUSION MATRIX
        function obj = setDiffusion(obj, D)
            obj.D = MatrixFunction(D, obj.Nx, obj.Ny);
        end
        
        % SETTING THE BOUNDARY CONDITIONS AS THE TRACE OF THE GIVEN
        % FUNCTION
        function obj = setBoundaryConditions(obj, f)
            obj.boundary_condition1 =@(X) f(X, obj.y_min);
            obj.boundary_condition2 =@(Y) f(obj.x_max, Y);
            obj.boundary_condition3 =@(X) f(X, obj.y_max);
            obj.boundary_condition4 =@(Y) f(obj.x_min, Y);
        end
        
        % SETTING THE SECOND MEMBER OCCURING IN THE PRESSURE EQUATION
        function obj = setSecondMemberDarcy(obj, f)
            obj.f1 = f;
        end
        
        % SETTING THE SECOND MEMBER OCCURING IN THE CONVECTION EQUATION
        function obj = setSecondMemberConvection(obj, f)
            obj.f2 = f;
        end
        
        % SETTING THE INITIAL CONDITION
        function obj = setInitialCondition(obj, f)
            obj.initial_condition = f;
        end
        
        % ADD A SUBDOMAIN BY GIVING ITS GEMEOMETRY AND THE ASSOCIATED
        % VALUES OF THE PERMEABILITY AND DIFFUSION MATRIX
        function obj = addSubdomain(obj, domain, K, D)
            nb_subdomains = length(obj.subdomains);
            obj.subdomains{nb_subdomains+1} = domain;
            obj.K.addSubdomain(domain, K);
            obj.D.addSubdomain(domain, D);
        end

        % PLOT THE SUBDOMAIN GEOMETRY AT POSITION Z_POS ON THE Z-AXIS
        function showSubdomains(obj, color, z_pos)
            nb_subdomains = length(obj.subdomains);
            for i = 1 : nb_subdomains
                obj.subdomains{i}.show(color, z_pos);
            end
        end
        
        % GENERATE A MESH
        function [X, Y] = generateMesh(obj)
            dx = (obj.x_max - obj.x_min) / (obj.Nx-1);
            dy = (obj.y_max - obj.y_min) / (obj.Ny-1);
            [X, Y] = meshgrid(obj.x_min : dx : obj.x_max, ...
                              obj.y_min : dy : obj.y_max);
        end
    end

end

