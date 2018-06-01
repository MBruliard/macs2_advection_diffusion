classdef MatrixFunction < handle
    % MatrixFunction: creates a 2 by 2 matrix function with the possibility
    % to specify sub-domains associated to different values
    % (such as obstacles)

    %% PROPERTIES DECLARATION

    properties
        m_default_value      % Main domain 2 by 2 matrix
        m_default_value_     % Precomputed coefficients of the main domain
                             % matrix multiplied by the
                             % main_domain_geometry matrix (see below)
        m_subdomain_value    % 3D matrix (last dimension refers to the
                             % subdomain id, the remaining domensions are
                             % associated to a 2 by 2 matrix
        nb_subdomains        % number of subdomains
        subdomain_geometry   % 3D matrix (last dimension refers to the
                             % subdomain id, which is associated to a 2D
                             % matrix with value 1 if cell is inside the
                             % subdomain, 0 otherwise)
        main_domain_geometry % 2D matrix of 0/1 values (1 if cell is inside
                             % the main domain, 0 otherwise)
    end


    %% METHODS DEFINITION

    methods
        % CONSTRUCTOR
        function obj = MatrixFunction(varargin)
            obj.nb_subdomains = 0;
            obj.m_default_value = zeros(2,2);
            if 3 == nargin
                M = varargin{1};
                obj.m_default_value = M;
            elseif 6 == nargin
                for i = 1 : 2
                    for j = 1 : 2
                        obj.m_default_value(i,j) = varargin{2*(i-1)+j};
                    end
                end
            else
                error('Wrong number of argument to create instance of "MatrixFunction" object')
            end
            obj.main_domain_geometry = ones(varargin{nargin}, varargin{nargin-1});
            for i = 1:2
                for j = 1:2
                    obj.m_default_value_(:,:,i,j) = M(i,j) ...
                        * obj.main_domain_geometry;
                end
            end
        end

        % ABILITY TO SPECIFY A SUBDOMAIN AND THE ASSOCIATED FUNCTION
        function obj = addSubdomain(obj, domain, M)
            % We convert M to a cell array of "function_handle objects" if
            % necessary
            obj.nb_subdomains = obj.nb_subdomains + 1;
            obj.subdomain_geometry(:,:,obj.nb_subdomains) = domain.pos_matrix;
            obj.main_domain_geometry = obj.main_domain_geometry ...
                - obj.main_domain_geometry .* domain.pos_matrix;
            obj.m_subdomain_value(:,:,obj.nb_subdomains) = M;
            
            for i = 1:2
                for j = 1:2
                    obj.m_default_value_(:,:,i,j) = obj.m_default_value(i,j) ...
                        * obj.main_domain_geometry;
                end
            end

        end

        % EVALUATE THE (i,j) COEFFICIENT AT THE LEVEL OF INDEX POSITION
        % (x_,y_) WITH x_ AND/OR y_ VECTOR(S)
        function res = eval_coef(obj, i, j, x_, y_)
            for k = obj.nb_subdomains : -1 : 1
                if obj.subdomain_geometry(y_,x_,k)
                    res = obj.m_subdomain_value(i,j,k);
                    return
                end
            end
            res = obj.m_default_value(i,j);
            return
        end

        % [VECTORIZED VERSION OF eval_coef]
        % EVALUATE THE (i,j) COEFFICIENT AT THE LEVEL OF INDEX POSITION
        % (x_,y_) WITH x_ AND/OR y_ VECTOR(S)
        function res = eval_coef_(obj, i, j, x_, y_)
            res = obj.m_default_value_(y_,x_,i,j);

            for k = obj.nb_subdomains : -1 : 1
                res = res + obj.m_subdomain_value(i,j,k) * obj.subdomain_geometry(y_,x_,k);
            end
            return

        end

        % EVALUATE THE HARMONIC MEAN OF THE (i,j) COEFFICIENT AT THE LEVEL
        % OF INTERFACE (x_index,y_index) + delta/2
        function res = eval_harmonic_mean(obj, i, j, x_index, y_index, delta)
            current_value  = obj.eval_coef(i, j, x_index, y_index);
            neighbor_value = obj.eval_coef(i, j, x_index+delta(1), y_index+delta(2));
            if 0 == current_value || 0 == neighbor_value
                res = 0;
                return
            end
            res = 2 * current_value * neighbor_value ...
                  / (current_value + neighbor_value);
        end

        % [VECTORIZED VERSION OF eval_harmonic_mean]
        % EVALUATE THE HARMONIC MEAN OF THE (i,j) COEFFICIENT AT THE LEVEL
        % OF INTERFACE (x_index,y_index) + delta/2
        function res = eval_harmonic_mean_(obj, i, j, x_, y_, delta)
            current_value  = obj.eval_coef_(i, j, x_, y_);
            neighbor_value = obj.eval_coef_(i, j, x_+delta(1), y_+delta(2));
            res = 2 * current_value .* neighbor_value ...
                  ./ (current_value + neighbor_value + ((0 == current_value & 0 == neighbor_value))) ...
                  .* (1 - (0 == current_value & 0 == neighbor_value));
        end

        % GET THE MAXIMUM ABSOLUTE COEFFICIENT AT POSITION (x,y)
        function m = get_max_coef_at_position(obj, x_index, y_index)

            for k = obj.nb_subdomains : -1 : 1
                if 1 == obj.subdomain_geometry(y_index, x_index, k)
                    m = max(max(abs(obj.m_subdomain_value(:,:,k))));
                    return
                end
            end
 
            m = max(max(abs(obj.m_default_value(:,:))));
        end

        % GET THE MAXIMUM ABSOLUTE COEFFICIENT
        function m = get_max_coef(obj, x_index, y_index)
            m = obj.get_max_coef_at_position(x_index(1), y_index(1));
            for i = x_index
                for j = y_index
                    l = obj.get_max_coef_at_position(i, j);
                    if l > m
                        m = l;
                    end
                end
            end
        end

    end

end

