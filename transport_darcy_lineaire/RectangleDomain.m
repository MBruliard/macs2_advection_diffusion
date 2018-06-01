classdef RectangleDomain
    % RectangleDomain: defines a rectangular domain object
    % A RectangleDomain object is caracterized by its lower left and
    % upper right points coordinates

    %% PROPERTIES DECLARATION

    properties
        pointA
        pointB
        pos_matrix
    end


    %% METHODS DEFINITION

    methods
        % CONSTRUCTOR
        function obj = RectangleDomain(varargin)
            if 4 == nargin
                obj.pointA.x = varargin{1}.x;
                obj.pointA.y = varargin{1}.y;
                obj.pointB.x = varargin{2}.x;
                obj.pointB.y = varargin{2}.y;
                last_arg = 3;
            elseif 6 == nargin
                obj.pointA.x = varargin{1};
                obj.pointA.y = varargin{2};
                obj.pointB.x = varargin{3};
                obj.pointB.y = varargin{4};
                last_arg = 5;
            else
                error('Wrong number of argument to create instance of "RectangleDomain" object')
            end
            
            GridX = varargin{last_arg};
            GridY = varargin{last_arg+1};
            Nx = size(GridX, 2);
            Ny = size(GridX, 1);
            obj.pos_matrix = zeros(Ny, Nx);
            for j = 1:Ny
                for i = 1:Nx
                    if    obj.pointA.x <= GridX(j,i) && GridX(j,i) <= obj.pointB.x ...
                       && obj.pointA.y <= GridY(j,i) && GridY(j,i) <= obj.pointB.y
                        obj.pos_matrix(j,i) = 1;
                    end
                end
            end
            
        end
        
        % CHECK IF A POINT IS CONTAINED IN THE DOMAIN GIVEN ITS COORDINATES
        function bool = contains_pos(obj, i, j)
            bool = obj.pos_matrix(j,i);
        end

        % CHECK IF A POINT IS CONTAINED IN THE DOMAIN GIVEN ITS COORDINATES
        function bool = contains_coord(obj, x, y)
            if    obj.pointA.x <= x && x <= obj.pointB.x ...
               && obj.pointA.y <= y && y <= obj.pointB.y
                bool = true;
            else
                bool = false;
            end
        end

        % CHECK IF A POINT IS CONTAINED IN THE DOMAIN GIVEN ITS VALUE
        function bool = contains_point(obj, point)
            if    obj.pointA.x <= point.x && point.x <= obj.pointB.x ...
               && obj.pointA.y <= point.y && point.y <= obj.pointB.y
                bool = true;
            else
                bool = false;
            end
        end

        % PLOT THE RECTANGULAR DOMAIN WITH THE INDICATED COLOR
        function show(obj, color, z_pos)
            patch([obj.pointA.x, obj.pointB.x, obj.pointB.x, obj.pointA.x, obj.pointA.x], ...
                     [obj.pointA.y, obj.pointA.y, obj.pointB.y, obj.pointB.y, obj.pointA.y], ...
                     [z_pos, z_pos, z_pos, z_pos, z_pos], color);
        end
    end

end
