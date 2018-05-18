classdef MatrixFunction < handle
    % MatrixFunction: creates a 2 by 2 matrix function with the possibility
    % to specify sub-domains associated to different values
    % (such as obstacles)

    %% PROPERTIES DECLARATION

    properties
        m_default_value
        m_subdomain_value
        subdomains
    end


    %% METHODS DEFINITION

    methods
        % Constructor
        function obj = MatrixFunction(varargin)
            obj.m_subdomain_value = {};
            obj.subdomains= {};
            if 1 == nargin
                % If the argument is an array we convert it into a cell array
                if isa(varargin{1}, 'cell')
                    M = varargin{1};
                elseif ismatrix(varargin{1})
                    M = num2cell(varargin{1});
                else
                    error('Wrong type given to MatrixFunction constructor')
                end
                % We convert all constant values into "function_handle"
                % objects
                N = cell(2,2);
                for i = 1 : 2
                    for j = 1 : 2
                        if isa(M{i,j}, 'function_handle')
                            N{i,j} = M{i,j};
                        else
                            N{i,j} =@(x,y) M{i,j};
                        end
                    end
                end
                obj.m_default_value = N;
            elseif 4 == nargin
                % Four arguments given, each corresponding to a coefficient
                % function
                for i = 1 : 2
                    for j = 1 : 2
                        if isa(varargin{2*(i-1)+j}, 'finction_handle')
                            obj.m_default_value{i,j} = varargin{2*(i-1)+j};
                        else
                            obj.m_default_value{i,j} =@(x,y) varargin{2*(i-1)+j};
                        end
                    end
                end
            else
                error('Wrong number of argument to create instance of "MatrixFunction" object')
            end
        end

        % Ability to specify a sub-domain and the associated function
        function obj = addSubdomain(obj, domain, M)
            % We convert M to a cell array of "function_handle objects" if
            % necessary
            if ~isa(M, 'cell')
                M = num2cell(M);
            end
            N = cell(2,2);
            for i = 1 : 2
                for j = 1 : 2
                    if isa(M{i,j}, 'function_handle')
                        N{i,j} = M{i,j};
                    else
                        N{i,j} =@(x,y) M{i,j};
                    end
                end
            end
            nb_subdomains = length(obj.subdomains);
            obj.subdomains{nb_subdomains+1} = domain;
            obj.m_subdomain_value{nb_subdomains+1} = N;
        end

        % Evaluate the (i,j) coefficient at the level of position (x,y)
        function res = eval_coef(obj, i, j, x, y)
            nb_subdomains = length(obj.subdomains);
            for k = nb_subdomains : -1 : 1
                if obj.subdomains{k}.contains_coord(x,y)
                    res = obj.m_subdomain_value{k}{i,j}(x,y);
                    return
                end
            end
            res = obj.m_default_value{i,j}(x,y);
        end

        % Evaluate the whole matrix at the level of position (x,y)
        function M = eval_matrix(obj, x, y)
            nb_subdomains = length(obj.subdomains);
            for k = nb_subdomains : -1 : 1
                if obj.subdomains{k}.contains_coord(x,y)
                    M(1,1) = obj.m_obstacle_value{k}{1,1}(x,y);
                    M(1,2) = obj.m_obstacle_value{k}{1,2}(x,y);
                    M(2,1) = obj.m_obstacle_value{k}{2,1}(x,y);
                    M(2,2) = obj.m_obstacle_value{k}{2,2}(x,y);
                    return
                end
            end
            M(1,1) = obj.m_default_value{1,1}(x,y);
            M(1,2) = obj.m_default_value{1,2}(x,y);
            M(2,1) = obj.m_default_value{2,1}(x,y);
            M(2,2) = obj.m_default_value{2,2}(x,y);
        end
    end

end

